"""A work queue on a shared filesystem: PBS jobs plus backfill from any machine.

    python validation/workqueue.py init    DIR TASKS.jsonl
    python validation/workqueue.py work    DIR [--workers N] [--hours H] [--nice] [--stale MIN]
    python validation/workqueue.py status  DIR
    python validation/workqueue.py collect DIR OUT.jsonl

TASKS.jsonl has one task per line, {"fn": "module:function", "args": {...}};
fn is imported from validation/ (validation/tasks.py has the usual ones) and
called as fn(**args), returning a JSON-serialisable result.

Any number of workers on any machines that see DIR (PBS jobs on headnode,
orr, cartman, stan: all mount /import/taiji1) pull tasks until none are left:

  * a task is claimed by creating DIR/claims/<i> (mkdir is atomic, NFS too);
  * its result goes to DIR/done/<i>.json, written to a temporary file and
    renamed, so a result is either complete or absent;
  * a worker touches its claim every minute while the task runs; a claim
    with no result whose heartbeat is older than --stale minutes (default
    10) is taken over, so tasks of workers killed by walltime, a dropped
    connection or a reboot are rerun. Tasks are deterministic, so running
    one twice is harmless.

A worker stops when every task is done or being worked on by a live worker
(heartbeat less than HEARTBEAT_MISSED old); while some unfinished task's
heartbeat is missed it waits, and takes the task over once its claim is
stale. With --hours it starts no task after that many hours (walltime).
"""

import argparse
import importlib
import json
import multiprocessing
import os
import socket
import sys
import threading
import time
import traceback

HEARTBEAT = 60            # seconds between claim touches
HEARTBEAT_MISSED = 150    # a claim this old belongs to a dead worker (probably)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path[:0] = [ROOT, HERE]


def _paths(qdir):
    return (os.path.join(qdir, 'tasks.jsonl'), os.path.join(qdir, 'claims'),
            os.path.join(qdir, 'done'))


def init(qdir, tasks_path):
    tasks_file, claims, done = _paths(qdir)
    if os.path.exists(tasks_file):
        sys.exit('%s already exists' % tasks_file)
    tasks = [json.loads(line) for line in open(tasks_path) if line.strip()]
    for t in tasks:
        assert ':' in t['fn'] and isinstance(t.get('args', {}), dict), t
    os.makedirs(claims)
    os.makedirs(done)
    with open(tasks_file, 'w') as f:
        for t in tasks:
            f.write(json.dumps(t) + '\n')
    print('queue %s: %d tasks' % (qdir, len(tasks)))


def _claim(claims, done, i, stale):
    c = os.path.join(claims, str(i))
    try:
        os.mkdir(c)
        return c
    except FileExistsError:
        pass
    try:
        age = time.time() - os.stat(c).st_mtime
    except FileNotFoundError:
        return None
    if age < stale or os.path.exists(os.path.join(done, '%d.json' % i)):
        return None
    # take over a stale claim: only one renamer can win
    try:
        os.rename(c, '%s.stale.%s.%d.%d' % (c, socket.gethostname(), os.getpid(), time.time()))
        os.mkdir(c)
        return c
    except OSError:
        return None


def _heartbeat(path, stop):
    while not stop.wait(HEARTBEAT):
        try:
            os.utime(path)
        except OSError:
            pass


def _run(task):
    mod, _, name = task['fn'].partition(':')
    fn = getattr(importlib.import_module(mod), name)
    return fn(**task.get('args', {}))


def _worker(qdir, deadline, stale, nice, seed):
    if nice:
        os.nice(19)
    tasks_file, claims, done = _paths(qdir)
    tasks = [json.loads(line) for line in open(tasks_file)]
    host = socket.gethostname().split('.')[0]
    n = len(tasks)
    # start at different offsets so workers rarely contend for the same claim
    order = [(seed * 7919 + k) % n for k in range(n)] if n else []
    while True:
        claimed = False
        for i in order:
            if deadline and time.time() > deadline:
                return
            out = os.path.join(done, '%d.json' % i)
            if os.path.exists(out):
                continue
            c = _claim(claims, done, i, stale)
            if c is None:
                continue
            claimed = True
            stop = threading.Event()
            threading.Thread(target=_heartbeat, args=(c, stop), daemon=True).start()
            t0 = time.time()
            rec = dict(task=i, host=host, pid=os.getpid())
            try:
                rec['result'] = _run(tasks[i])
            except Exception as e:   # record, don't abort the queue
                rec['error'] = '%s: %s' % (type(e).__name__, e)
                rec['traceback'] = traceback.format_exc()
            rec['seconds'] = round(time.time() - t0, 2)
            stop.set()
            tmp = '%s.tmp.%s.%d' % (out, host, os.getpid())
            with open(tmp, 'w') as f:
                f.write(json.dumps(rec) + '\n')
            os.replace(tmp, out)
        if claimed:
            continue
        # nothing claimable: wait only for tasks whose worker seems dead
        orphans = [i for i, age in _unfinished(claims, done) if age > HEARTBEAT_MISSED]
        if not orphans or (deadline and time.time() > deadline):
            return
        time.sleep(30)


def _unfinished(claims, done):
    """(task, heartbeat age in seconds) of claimed tasks without a result."""
    out, now = [], time.time()
    for d in os.listdir(claims):
        if d.isdigit() and not os.path.exists(os.path.join(done, d + '.json')):
            try:
                out.append((int(d), now - os.stat(os.path.join(claims, d)).st_mtime))
            except FileNotFoundError:
                pass
    return out


def work(qdir, workers, hours, stale, nice):
    deadline = time.time() + hours * 3600 if hours else None
    stale = max(stale, HEARTBEAT_MISSED / 60)   # never take over a live claim
    procs = [multiprocessing.Process(target=_worker, args=(qdir, deadline, stale * 60, nice, k + os.getpid()))
             for k in range(workers)]
    for p in procs:
        p.start()
    for p in procs:
        p.join()


def status(qdir):
    tasks_file, claims, done = _paths(qdir)
    n = sum(1 for _ in open(tasks_file))
    results = [json.loads(open(os.path.join(done, f)).read())
               for f in os.listdir(done) if f.endswith('.json')]
    active = _unfinished(claims, done)
    running = [i for i, age in active if age <= HEARTBEAT_MISSED]
    hosts = {}
    for r in results:
        hosts[r['host']] = hosts.get(r['host'], 0) + 1
    errors = sum('error' in r for r in results)
    cpu = sum(r['seconds'] for r in results)
    print('%d tasks: %d done (%d errors), %d running, %d orphaned, %d to go; %.0f CPU-s so far'
          % (n, len(results), errors, len(running), len(active) - len(running),
             n - len(results) - len(running), cpu))
    print('done by host:', ', '.join('%s %d' % kv for kv in sorted(hosts.items())))


def collect(qdir, out_path):
    tasks_file, _, done = _paths(qdir)
    n = sum(1 for _ in open(tasks_file))
    missing = 0
    with open(out_path, 'w') as f:
        for i in range(n):
            p = os.path.join(done, '%d.json' % i)
            if not os.path.exists(p):
                missing += 1
                continue
            f.write(open(p).read())
    print('%d results written, %d missing' % (n - missing, missing))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', choices=('init', 'work', 'status', 'collect'))
    ap.add_argument('dir')
    ap.add_argument('path', nargs='?')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--hours', type=float, default=None,
                    help='do not start tasks after this many hours (PBS walltime margin)')
    ap.add_argument('--stale', type=float, default=10, help='minutes without heartbeat')
    ap.add_argument('--nice', action='store_true', help='run at the lowest priority (backfill)')
    a = ap.parse_args()
    if a.cmd == 'init':
        init(a.dir, a.path)
    elif a.cmd == 'work':
        work(a.dir, a.workers, a.hours, a.stale, a.nice)
    elif a.cmd == 'status':
        status(a.dir)
    else:
        collect(a.dir, a.path)


if __name__ == '__main__':
    main()
