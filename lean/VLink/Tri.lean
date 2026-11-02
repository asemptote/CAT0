import Mathlib

/-!
# Closed triangulated 3-manifolds from isomorphism signatures (layer L3)

Mirrors `vlink/tri.py`: decodes Regina's isomorphism signatures into a gluing
table, numbers edges and vertices in order of first appearance over
(tetrahedron, local index) as Regina does, and rejects anything that is not a
closed 3-manifold triangulation (boundary facets, an edge identified with
itself in reverse, a vertex link that is not a 2-sphere).

`adj[t][f]` is the tetrahedron glued to face `f` of `t` (the face opposite
vertex `f`); `glue[t][f]` the vertex map: vertex `i` of `t` goes to vertex
`glue[t][f][i]` of `adj[t][f]`.
-/

namespace VLink

def b64 (c : Char) : Option Nat :=
  if 'a' ≤ c ∧ c ≤ 'z' then some (c.toNat - 'a'.toNat)
  else if 'A' ≤ c ∧ c ≤ 'Z' then some (c.toNat - 'A'.toNat + 26)
  else if '0' ≤ c ∧ c ≤ '9' then some (c.toNat - '0'.toNat + 52)
  else if c = '+' then some 62
  else if c = '-' then some 63
  else none

/-- Regina's `Perm<4>::orderedSn`: permutations of `{0,1,2,3}` in lexicographic order. -/
def perms4 : Array (Array Nat) := Id.run do
  let mut out := #[]
  for a in [0, 1, 2, 3] do
    for b in [0, 1, 2, 3] do
      for c in [0, 1, 2, 3] do
        if a ≠ b ∧ a ≠ c ∧ b ≠ c then out := out.push #[a, b, c, 6 - a - b - c]
  return out

def inv4 (p : Array Nat) : Array Nat := Id.run do
  let mut q := #[0, 0, 0, 0]
  for i in [0:4] do q := q.set! p[i]! i
  return q

/-- Regina's local edge numbering in a tetrahedron. -/
def EDGES : Array (Nat × Nat) := #[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

def edgeIndex (a b : Nat) : Nat :=
  (EDGES.findIdx? (fun e => e = (a, b) ∨ e = (b, a))).getD 0

/-- Read an `nchars`-character little-endian base-64 integer at `pos`. -/
def readInt (s : Array Char) (pos nchars : Nat) : Except String (Nat × Nat) := do
  let mut v := 0
  for i in [0:nchars] do
    match s[pos + i]? >>= b64 with
    | some x => v := v ||| (x <<< (6 * i))
    | none => throw "incomplete or invalid base64 encoding"
  return (v, pos + nchars)

/-- The gluing table of a closed triangulation (Regina's `fromIsoSig`). -/
def decode (sig : String) : Except String (Array (Array Nat) × Array (Array (Array Nat))) := do
  let s := sig.trimAscii.toString.toList.toArray
  let mut pos := 0
  let mut adj : Array (Array (Option Nat)) := #[]
  let mut glue : Array (Array (Array Nat)) := #[]
  while pos < s.size do
    let (n0, p1) ← readInt s pos 1
    pos := p1
    let mut n := n0
    let mut nchars := 1
    if n0 = 63 then
      let (nc, p2) ← readInt s pos 1
      let (nn, p3) ← readInt s p2 nc
      pos := p3
      nchars := nc
      n := nn
    if n = 0 then continue
    let mut actions : Array Nat := #[]
    let mut nf := 0
    let mut nj := 0
    while nf < 4 * n do
      let (c, p) ← readInt s pos 1
      pos := p
      for j in [0:3] do
        let a := (c >>> (2 * j)) &&& 3
        if nf = 4 * n then
          if a ≠ 0 then throw "extraneous facet actions"
        else
          if a = 3 then throw "invalid facet action"
          nf := nf + (if a = 0 then 1 else 2)
          if a = 2 then nj := nj + 1
          if nf > 4 * n then throw "facet actions do not match triangulation size"
          actions := actions.push a
    let mut dest := #[]
    for _ in [0:nj] do
      let (d, p) ← readInt s pos nchars
      pos := p
      dest := dest.push d
    let mut perm := #[]
    for _ in [0:nj] do
      let (d, p) ← readInt s pos 1
      pos := p
      perm := perm.push d
    if s[pos]? = some '.' then throw "lock masks are not supported"
    let base := adj.size
    for _ in [0:n] do
      adj := adj.push #[none, none, none, none]
      glue := glue.push #[#[], #[], #[], #[]]
    let mut k := 0
    let mut nxt := 1
    let mut jp := 0
    for p in [0:n] do
      for f in [0:4] do
        if (adj[base + p]!)[f]!.isSome then continue
        let a := actions[k]!
        k := k + 1
        if a = 0 then throw "boundary facet: the triangulation is not closed"
        let mut q := 0
        let mut g := #[0, 1, 2, 3]
        if a = 1 then
          if nxt ≥ n then throw "gluing to non-existent simplex"
          q := nxt
          nxt := nxt + 1
        else
          if perm[jp]! ≥ 24 then throw "invalid gluing permutation"
          q := dest[jp]!
          g := perms4[perm[jp]!]!
          jp := jp + 1
          if q ≥ nxt || (adj[base + q]!)[g[f]!]!.isSome || (q = p && g[f]! = f) then
            throw "invalid gluing destination"
        adj := adj.modify (base + p) (·.set! f (some (base + q)))
        glue := glue.modify (base + p) (·.set! f g)
        adj := adj.modify (base + q) (·.set! g[f]! (some (base + p)))
        glue := glue.modify (base + q) (·.set! g[f]! (inv4 g))
  if adj.size = 0 then throw "empty triangulation"
  return (adj.map (·.map (·.getD 0)), glue)

/-- A closed triangulated 3-manifold with its skeleton. -/
structure Tri where
  n : Nat
  adj : Array (Array Nat)
  glue : Array (Array (Array Nat))
  /-- `edgeOf[t][k]`: the edge of local edge `k` of `t`. -/
  edgeOf : Array (Array Nat)
  /-- `edges[e]`: the sorted `(t, k)` of edge `e`. -/
  edges : Array (Array (Nat × Nat))
  /-- `ends[(t, v, w)] = (e, 0 or 1)`: which end of its edge the vertex `v` of edge `{v, w}` is. -/
  ends : Std.HashMap (Nat × Nat × Nat) (Nat × Nat)
  vertexOf : Array (Array Nat)
  /-- `corners[i]`: the sorted corners `(t, v)` of vertex `i`. -/
  corners : Array (Array (Nat × Nat))
  deriving Inhabited

def sortPairs (l : Array (Nat × Nat)) : Array (Nat × Nat) :=
  l.qsort (fun a b => a.1 < b.1 ∨ (a.1 = b.1 ∧ a.2 < b.2))

def Tri.ofIsoSig (sig : String) : Except String Tri := do
  let (adj, glue) ← decode sig
  let n := adj.size
  -- edges
  let mut edgeOf : Array (Array (Option Nat)) := Array.replicate n (Array.replicate 6 none)
  let mut edges : Array (Array (Nat × Nat)) := #[]
  let mut ends : Std.HashMap (Nat × Nat × Nat) (Nat × Nat) := {}
  for t in [0:n] do
    for k in [0:6] do
      if (edgeOf[t]!)[k]!.isSome then continue
      let e := edges.size
      let (a0, b0) := EDGES[k]!
      let mut seen : Std.HashSet (Nat × Nat × Nat) := ({} : Std.HashSet _).insert (t, a0, b0)
      let mut stack := #[(t, a0, b0)]
      let mut members := #[]
      edgeOf := edgeOf.modify t (·.set! k (some e))
      while h : stack.size > 0 do
        let (u, a, b) := stack.back
        stack := stack.pop
        members := members.push (u, edgeIndex a b)
        ends := (ends.insert (u, a, b) (e, 0)).insert (u, b, a) (e, 1)
        for f in [0:4] do
          if f = a ∨ f = b then continue
          let g := (glue[u]!)[f]!
          let w := (adj[u]!)[f]!
          let (a', b') := (g[a]!, g[b]!)
          if seen.contains (w, b', a') then
            throw s!"edge {e} is identified with itself in reverse"
          if !seen.contains (w, a', b') then
            seen := seen.insert (w, a', b')
            edgeOf := edgeOf.modify w (·.set! (edgeIndex a' b') (some e))
            stack := stack.push (w, a', b')
      edges := edges.push (sortPairs members)
  -- vertices
  let mut vertexOf : Array (Array (Option Nat)) := Array.replicate n (Array.replicate 4 none)
  let mut corners : Array (Array (Nat × Nat)) := #[]
  for t in [0:n] do
    for v in [0:4] do
      if (vertexOf[t]!)[v]!.isSome then continue
      let i := corners.size
      vertexOf := vertexOf.modify t (·.set! v (some i))
      let mut stack := #[(t, v)]
      let mut members := #[]
      while h : stack.size > 0 do
        let (u, x) := stack.back
        stack := stack.pop
        members := members.push (u, x)
        for f in [0:4] do
          if f = x then continue
          let w := (adj[u]!)[f]!
          let y := ((glue[u]!)[f]!)[x]!
          if (vertexOf[w]!)[y]!.isNone then
            vertexOf := vertexOf.modify w (·.set! y (some i))
            stack := stack.push (w, y)
      corners := corners.push (sortPairs members)
  -- every vertex link is a 2-sphere: V - E + F = 2
  for i in [0:corners.size] do
    let cs := corners[i]!
    let mut lvs : Std.HashSet (Nat × Nat) := {}
    for (t, v) in cs do
      for w in [0:4] do
        if w ≠ v then lvs := lvs.insert (ends.getD (t, v, w) (0, 0))
    if lvs.size + cs.size ≠ 2 + 3 * cs.size / 2 then
      throw s!"the link of vertex {i} is not a 2-sphere"
  return {
    n, adj, glue,
    edgeOf := edgeOf.map (·.map (·.getD 0)),
    edges, ends,
    vertexOf := vertexOf.map (·.map (·.getD 0)),
    corners }

/-- Link-triangle vertex `j` ↦ tetrahedron vertex, for the corner at `v`:
the other vertices in increasing order, then `v`. -/
def cornerPerm (v : Nat) : Array Nat :=
  ((#[0, 1, 2, 3].filter (· ≠ v)).push v)

end VLink
