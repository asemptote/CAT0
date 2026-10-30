import VLink.Vec
import VLink.Tri

/-!
# The kernel's geometry and certificate checks (layers L3–L4, executable)

A transcription of `vlink/kernel.py` with every decision exact: signs by
`MQ.sign` / `Tw.sign` (proved sound), angle sums by `AngleSum` (proved
correct). Functions return `Except String`: an error is a rejection, either
because a claim is false or because a sign ran out of fuel.

Unlike the Python kernel, every length comparison is exact (Python's
`exact_only=True`).
-/

namespace VLink

def twSgn (x : Tw) : KM Int :=
  match Tw.sign kfuel x with
  | some s => pure s
  | none => throw "sign undecided (out of fuel)"

def mqOne : MQ := MQ.ofRat 1

/-! ### Vectors: `V3` and the predicates on it are in `VLink/Vec.lean` -/

abbrev Vec := V3

def E3 : Array V3 := #[⟨mqOne, [], []⟩, ⟨[], mqOne, []⟩, ⟨[], [], mqOne⟩]

/-! ### Angle sums -/

def angleAdd (st : AngleSum) (ab : MQ × MQ) : KM AngleSum :=
  match st.add ab.1 ab.2 with
  | some st' => pure st'
  | none => throw "angle undecided or between degenerate vectors"

def angleSumL : AngleSum → List (MQ × MQ) → KM AngleSum
  | st, [] => .ok st
  | st, p :: l =>
    match st.add p.1 p.2 with
    | some st' => angleSumL st' l
    | none => .error "angle undecided or between degenerate vectors"

def angleSum (l : Array (MQ × MQ)) : KM AngleSum := angleSumL AngleSum.zero l.toList

/-- The exact angle sum represents the sum of the angles. -/
theorem angleSumL_ok {st st' : AngleSum} {Θ : ℝ} {l : List (MQ × MQ)} (hr : st.Rep Θ)
    (h : angleSumL st l = .ok st') :
    st'.Rep (Θ + (l.map fun p => angleOf p.1 p.2).sum) := by
  induction l generalizing st Θ with
  | nil => cases h; simpa using hr
  | cons p l ih =>
    unfold angleSumL at h
    split at h
    · rename_i st1 h1
      have := ih (AngleSum.add_sound hr h1) h
      simpa [add_assoc] using this
    · cases h

theorem angleSum_ok {l : Array (MQ × MQ)} {st : AngleSum} (h : angleSum l = .ok st) :
    st.Rep ((l.toList.map fun p => angleOf p.1 p.2).sum) := by
  have := angleSumL_ok AngleSum.rep_zero h
  simpa using this

def atLeast (st : AngleSum) (k : Nat) : Bool := k ≤ st.h

def angleLeq (x y : AngleSum) : KM Bool :=
  match x.leq y with
  | some b => pure b
  | none => throw "angle comparison undecided"

/-! ### Metric -/

def gramDet (rows : Gram) : ℚ := Id.run do
  let mut M := rows
  let n := M.size
  let mut det : ℚ := 1
  for c in [0:n] do
    match (List.range' c (n - c)).find? (fun r => (M[r]!)[c]! ≠ 0) with
    | none => return 0
    | some p =>
      if p ≠ c then
        let rc := M[c]!
        M := (M.set! c M[p]!).set! p rc
        det := -det
      det := det * (M[c]!)[c]!
      for r in [c + 1:n] do
        let f := (M[r]!)[c]! / (M[c]!)[c]!
        let rowc := M[c]!
        M := M.modify r fun row => Id.run do
          let mut row := row
          for k in [c:n] do row := row.set! k (row[k]! - f * rowc[k]!)
          return row
  return det

structure Geometry where
  T : Tri
  /-- `d[t][i][j]`: squared length of the edge `ij` of tetrahedron `t`. -/
  d : Array (Array (Array ℚ))
  deriving Inhabited

def Geometry.make (T : Tri) (sqlens : Array ℚ) : KM Geometry := do
  if sqlens.size ≠ T.edges.size then throw "wrong number of edge lengths"
  let mut d := #[]
  for t in [0:T.n] do
    let mut m := Array.replicate 4 (Array.replicate 4 (0 : ℚ))
    for k in [0:6] do
      let (a, b) := EDGES[k]!
      let L := sqlens[(T.edgeOf[t]!)[k]!]!
      m := m.modify a (·.set! b L)
      m := m.modify b (·.set! a L)
    d := d.push m
  return ⟨T, d⟩

def Geometry.gram (g : Geometry) (tet base : Nat) (others : Array Nat) : Gram :=
  let d := g.d[tet]!
  others.map fun x => others.map fun y =>
    ((d[base]!)[x]! + (d[base]!)[y]! - (if x ≠ y then (d[x]!)[y]! else 0)) / 2

def Geometry.tetValid (g : Geometry) (tet : Nat) : Bool :=
  gramDet (g.gram tet 0 #[1, 2, 3]) > 0 &&
    [#[0, 1, 2], #[0, 1, 3], #[0, 2, 3], #[1, 2, 3]].all fun f =>
      gramDet (g.gram tet f[0]! #[f[1]!, f[2]!]) > 0

/-- `(a, b²)` of the dihedral angle at edge `(i, j)` of `tet`. -/
def Geometry.dihedral (g : Geometry) (tet i j : Nat) : MQ × MQ :=
  let kl := #[0, 1, 2, 3].filter (fun x => x ≠ i ∧ x ≠ j)
  let G := g.gram tet i #[j, kl[0]!, kl[1]!]
  let ee := (G[0]!)[0]!
  let eu := (G[0]!)[1]!
  let ew := (G[0]!)[2]!
  let a := (G[1]!)[2]! - eu * ew / ee
  let nu := (G[1]!)[1]! - eu * eu / ee
  let nw := (G[2]!)[2]! - ew * ew / ee
  (MQ.ofRat a, MQ.ofRat (nu * nw - a * a))

def Geometry.edgeAngle (g : Geometry) (e : Nat) : KM AngleSum :=
  angleSum (g.T.edges[e]!.map fun (tet, k) => g.dihedral tet EDGES[k]!.1 EDGES[k]!.2)

/-! ### Vertex links -/

structure Trans where
  t2 : Nat
  g : Array Nat
  p : ℚ × ℚ
  q : ℚ × ℚ
  mu : MQ
  deriving Inhabited

structure Link where
  n : Nat
  G : Array Gram
  adj : Array (Array Nat)
  glue : Array (Array (Array Nat))
  lv : Array (Array Nat)
  nVertices : Nat
  trans : Array (Array Trans)
  deriving Inhabited

def others3 (d : Nat) : Nat × Nat := if d = 0 then (1, 2) else if d = 1 then (0, 2) else (0, 1)

def mkTrans (G : Array Gram) (adj : Array (Array Nat)) (glue : Array (Array (Array Nat)))
    (t d : Nat) : KM Trans := do
  let t2 := (adj[t]!)[d]!
  let g := (glue[t]!)[d]!
  let (a, b) := others3 d
  let Gt := G[t]!
  let G2 := G[t2]!
  let (a2, b2, d2) := (g[a]!, g[b]!, g[d]!)
  if (Gt[a]!)[a]! ≠ (G2[a2]!)[a2]! ∨ (Gt[a]!)[b]! ≠ (G2[a2]!)[b2]! ∨ (Gt[b]!)[b]! ≠ (G2[b2]!)[b2]! then
    throw "glued link edge has different lengths on its two sides"
  let dt := (Gt[a]!)[a]! * (Gt[b]!)[b]! - (Gt[a]!)[b]! ^ 2
  let p := (((Gt[d]!)[a]! * (Gt[b]!)[b]! - (Gt[d]!)[b]! * (Gt[a]!)[b]!) / dt,
            ((Gt[d]!)[b]! * (Gt[a]!)[a]! - (Gt[d]!)[a]! * (Gt[a]!)[b]!) / dt)
  let q := (((G2[d2]!)[a2]! * (Gt[b]!)[b]! - (G2[d2]!)[b2]! * (Gt[a]!)[b]!) / dt,
            ((G2[d2]!)[b2]! * (Gt[a]!)[a]! - (G2[d2]!)[a2]! * (Gt[a]!)[b]!) / dt)
  let h := (Gt[d]!)[d]! - p.1 * (Gt[d]!)[a]! - p.2 * (Gt[d]!)[b]!
  let h2 := (G2[d2]!)[d2]! - q.1 * (G2[d2]!)[a2]! - q.2 * (G2[d2]!)[b2]!
  if h ≤ 0 ∨ h2 ≤ 0 then throw "degenerate link triangle"
  return ⟨t2, g, p, q, MQ.sqrtRat (h2 / h)⟩

def Link.make (geo : Geometry) (iv : Nat) : KM Link := do
  let T := geo.T
  let corners := T.corners[iv]!
  let mut index : Std.HashMap (Nat × Nat) Nat := {}
  for i in [0:corners.size] do index := index.insert corners[i]! i
  let mut G := #[]
  let mut adj := #[]
  let mut glue := #[]
  let mut lv := #[]
  let mut labels : Std.HashMap (Nat × Nat) Nat := {}
  for (tet, v) in corners do
    let fp := cornerPerm v
    G := G.push (geo.gram tet v #[fp[0]!, fp[1]!, fp[2]!])
    let mut a := #[]
    let mut gl := #[]
    for d in [0:3] do
      let g := (T.glue[tet]!)[fp[d]!]!
      let v2 := g[v]!
      let fp2 := cornerPerm v2
      a := a.push (index.getD ((T.adj[tet]!)[fp[d]!]!, v2) 0)
      gl := gl.push (#[0, 1, 2].map fun j => (fp2.findIdx? (· = g[fp[j]!]!)).getD 0)
    adj := adj.push a
    glue := glue.push gl
    let mut l := #[]
    for j in [0:3] do
      let key := T.ends.getD (tet, v, fp[j]!) (0, 0)
      match labels.get? key with
      | some x => l := l.push x
      | none =>
        l := l.push labels.size
        labels := labels.insert key labels.size
    lv := lv.push l
  let mut trans := #[]
  for t in [0:corners.size] do
    let mut row := #[]
    for d in [0:3] do row := row.push (← mkTrans G adj glue t d)
    trans := trans.push row
  return ⟨corners.size, G, adj, glue, lv, labels.size, trans⟩

/-- Develop across the edge opposite vertex `d` of triangle `t` at positions `P`. -/
def step (L : Link) (P : Array Vec) (t d : Nat) : Nat × Array Vec :=
  let tr := (L.trans[t]!)[d]!
  let (a, b) := others3 d
  let Pa := P[a]!
  let Pb := P[b]!
  let Pd := P[d]!
  let normal := V3.sub (V3.sub Pd (V3.scaleQ tr.p.1 Pa)) (V3.scaleQ tr.p.2 Pb)
  let new := V3.sub (V3.add (V3.scaleQ tr.q.1 Pa) (V3.scaleQ tr.q.2 Pb)) (V3.scale tr.mu normal)
  let P2 := ((#[default, default, default] : Array V3).set! tr.g[a]! Pa |>.set! tr.g[b]! Pb).set!
    tr.g[d]! new
  (tr.t2, P2)

def develop (L : Link) (ts steps : Array Nat) : KM (Array (Array Vec)) := do
  if steps.size + 1 ≠ ts.size then throw "wrong number of strip steps"
  let mut P := E3
  let mut out := #[E3]
  for k in [0:steps.size] do
    let t := ts[k]!
    let d := steps[k]!
    if d > 2 ∨ (L.adj[t]!)[d]! ≠ ts[k + 1]! then
      throw s!"strip step {k} does not glue triangle {t} to {ts[k + 1]!}"
    P := (step L P t d).2
    out := out.push P
  return out

/-! ### Geodesics and angles -/

/-! `coneSigns`, `crosses`, `tangentAngle`, `length`: `VLink/Vec.lean` -/

/-! ### Beads -/

structure Bead where
  ts : Array Nat
  start : Nat
  stop : Nat
  G : Gram
  pos : Array (Array Vec)
  steps : Array Nat
  S : Vec
  E : Vec
  pi : Bool
  inDir : Vec
  outDir : Vec
  deriving Inhabited

def Bead.make (ts : Array Nat) (start stop : Nat) (G : Gram) (pos : Array (Array Vec))
    (steps : Array Nat) (Q : Option Vec) : KM Bead :=
  let S := (pos[0]!)[start]!
  let E := pos.back![stop]!
  match parallel G S E with
  | .error e => .error e
  | .ok par =>
    let pi := decide (par < 0)
    let q := Q.getD default
    .ok ⟨ts, start, stop, G, pos, steps, S, E, pi, if pi then q else S, if pi then q else E⟩

lemma Bead.make_ok {ts : Array Nat} {start stop : Nat} {G : Gram} {pos : Array (Array Vec)}
    {steps : Array Nat} {Q : Option Vec} {b : Bead}
    (h : Bead.make ts start stop G pos steps Q = .ok b) :
    b.ts = ts ∧ b.G = G ∧ b.pos = pos ∧ b.steps = steps ∧ b.start = start ∧ b.stop = stop ∧
      b.S = (pos[0]!)[start]! ∧ b.E = pos.back![stop]! ∧
      b.outDir = (if b.pi then Q.getD default else b.E) := by
  unfold Bead.make at h
  dsimp only at h
  split at h
  · cases h
  · cases h
    simp

/-! ### Bead enumeration (completeness: docs/proofs.md §3) -/

def visible (S A B : Vec) : KM (Option (Vec × Vec)) := do
  let o ← sgn (det3 S A B)
  if o == 0 then return none
  return some (if o > 0 then (A, B) else (B, A))

def inside (S L R X : Vec) (strict : Bool) : KM Bool := do
  let a ← sgn (det3 S L X)
  let b ← sgn (det3 S X R)
  return if strict then a > 0 && b > 0 else a ≥ 0 && b ≥ 0

def intersect (S : Vec) (I1 I2 : Vec × Vec) : KM (Option (Vec × Vec)) := do
  let (L1, R1) := I1
  let (L2, R2) := I2
  let L ← if ← inside S L1 R1 L2 false then pure (some L2)
    else if ← inside S L2 R2 L1 false then pure (some L1) else pure none
  let R ← if ← inside S L1 R1 R2 false then pure (some R2)
    else if ← inside S L2 R2 R1 false then pure (some R1) else pure none
  match L, R with
  | some L, some R => return if (← sgn (det3 S L R)) > 0 then some (L, R) else none
  | _, _ => return none

/-- Cross the edge opposite vertex `d` of the strip's last triangle (`_grow`). -/
def grow (L : Link) (G : Gram) (piBeads : Bool) :
    Nat → Array Nat → Array (Array Vec) → Array Nat → Nat → Vec → Nat → Vec × Vec →
      Array Bead → KM (Array Bead)
  | 0, _, _, _, _, _, _, _, _ => throw "bead strip longer than 4x the link: growth not terminating?"
  | fuel + 1, ts, pos, steps, start, S, d, iv, out => do
    let t := ts.back!
    let t2 := (L.adj[t]!)[d]!
    let (_, P2) := step L pos.back! t d
    let ts2 := ts.push t2
    let pos2 := pos.push P2
    let steps2 := steps.push d
    let r := ((L.glue[t]!)[d]!)[d]!
    let E := P2[r]!
    let par ← parallel G S E
    let mut out := out
    if par < 0 then
      if piBeads then out := out.push (← Bead.make ts2 start r G pos2 steps2 none)
    else if par == 0 then
      if ← inside S iv.1 iv.2 E true then
        out := out.push (← Bead.make ts2 start r G pos2 steps2 none)
    for dn in [(r + 1) % 3, (r + 2) % 3] do
      let y := 3 - r - dn
      if par != 0 then continue
      if (← parallel G S P2[y]!) != 0 then continue
      match ← visible S E P2[y]! with
      | none => continue
      | some iv2 =>
        match ← intersect S iv iv2 with
        | none => continue
        | some new => out ← grow L G piBeads fuel ts2 pos2 steps2 start S dn new out
    return out

def enumerateBeads (L : Link) (piBeads : Bool) : KM (Array Bead) := do
  let mut out := #[]
  for t in [0:L.n] do
    let G := L.G[t]!
    for d0 in [0:3] do
      let P := E3
      let a := (d0 + 1) % 3
      let b := (d0 + 2) % 3
      let S := P[d0]!
      out := out.push (← Bead.make #[t] d0 a G #[P] #[] none)
      out := out.push (← Bead.make #[t] d0 b G #[P] #[] none)
      if (← parallel G S P[a]!) != 0 || (← parallel G S P[b]!) != 0 then continue
      match ← visible S P[a]! P[b]! with
      | none => continue
      | some iv => out ← grow L G piBeads (4 * L.n + 1) #[t] #[P] #[] d0 S d0 iv out
  return out

/-! ### Junctions -/

/-- The angle at the junction `bin → bout` on side `k ∈ {1, 2}`. -/
def junctionSide (L : Link) (bin bout : Bead) (k : Nat) : KM AngleSum := do
  let tm := bin.ts.back!
  let rm := bin.stop
  let t0 := bout.ts[0]!
  let d0 := bout.start
  let P := bin.pos.back!
  let Gin := bin.G
  let W := P[rm]!
  let u := bin.inDir
  let d := (rm + k) % 3
  let y := 3 - rm - d
  if (tm, rm) == (t0, d0) then
    -- is the outgoing direction inside this corner, between u and ray(y)?
    let X : Vec := V3.add (V3.add (V3.scale bout.outDir.x P[0]!) (V3.scale bout.outDir.y P[1]!))
      (V3.scale bout.outDir.z P[2]!)
    let a1b1 := tangentAngle Gin W u X
    let s12 ← angleSum #[a1b1, tangentAngle Gin W X P[y]!]
    let s3 ← angleSum #[tangentAngle Gin W u P[y]!]
    if s12.h == s3.h then
      let im ← twSgn (Tw.sub (Tw.mul s12.im s3.re) (Tw.mul s12.re s3.im))
      if im == 0 then
        let re ← twSgn (Tw.add (Tw.mul s12.re s3.re) (Tw.mul s12.im s3.im))
        if re > 0 then return ← angleSum #[a1b1]
  let mut acc ← angleSum #[tangentAngle Gin W u P[y]!]
  let mut t := tm
  let mut c := rm
  let mut ent := d
  let mut ex := y
  for _ in [0:4 * L.n + 2] do
    if atLeast acc 1 then return acc
    let g := (L.glue[t]!)[ent]!
    t := (L.adj[t]!)[ent]!
    c := g[c]!
    ent := g[ex]!
    ex := 3 - c - ent
    if (t, c) == (t0, d0) then
      return ← angleAdd acc (tangentAngle bout.G E3[d0]! E3[ent]! bout.outDir)
    acc ← angleAdd acc (tangentAngle L.G[t]! E3[c]! E3[ent]! E3[ex]!)
  throw "junction fan walk did not reach the outgoing triangle"

/-- Angle `≥ π` on both sides of the junction `bin → bout`. -/
def junctionOk (L : Link) (bin bout : Bead) : KM Bool := do
  return atLeast (← junctionSide L bin bout 1) 1 && atLeast (← junctionSide L bin bout 2) 1

def beadLength (b : Bead) : MQ × MQ :=
  if b.pi then (MQ.ofRat (-1), []) else length b.G b.S b.E

def totalLength (beads : Array Bead) : KM AngleSum := angleSum (beads.map beadLength)

/-! ### Certificates -/

/-- A bead as given in a certificate: link triangles, start / end corners,
gluing steps, and for a π-bead a point fixing its geodesic. -/
structure BeadSpec where
  ts : Array Nat
  start : Nat
  stop : Nat
  steps : Array Nat
  piPoint : Option Vec
  deriving Inhabited

/-- Every interior edge `k` of the strip is crossed by the bead's geodesic. -/
def checkCrossings (G : Gram) (S Q : V3) (pi : Bool) (pos : Array (Array V3))
    (steps : Array Nat) : List Nat → KM Unit
  | [] => .ok ()
  | k :: ks =>
    match crosses G S Q ((pos[k]!)[(others3 steps[k]!).1]!) ((pos[k]!)[(others3 steps[k]!).2]!) pi with
    | .ok true => checkCrossings G S Q pi pos steps ks
    | .ok false => .error s!"geodesic does not cross interior edge {k}"
    | .error e => .error e

lemma checkCrossings_ok {G : Gram} {S Q : V3} {pi : Bool} {pos : Array (Array V3)}
    {steps : Array Nat} {l : List Nat} (h : checkCrossings G S Q pi pos steps l = .ok ()) :
    ∀ k ∈ l, crosses G S Q ((pos[k]!)[(others3 steps[k]!).1]!)
      ((pos[k]!)[(others3 steps[k]!).2]!) pi = .ok true := by
  induction l with
  | nil => simp
  | cons k l ih =>
    unfold checkCrossings at h
    split at h
    · rename_i hk
      intro j hj
      rcases List.mem_cons.mp hj with rfl | hj
      · exact hk
      · exact ih h j hj
    · cases h
    · cases h

/-- With a π-bead's point `q`: it must not be on the axis through `S`. -/
def piQ (G : Gram) (S : V3) : Option V3 → KM (Option V3)
  | none => .error "pi-bead without pi_point"
  | some q =>
    match parallel G S q with
    | .ok 0 => .ok (some q)
    | .ok _ => .error "pi_point on the axis"
    | .error e => .error e

/-- A successful step of a chain of checks. -/
lemma bind_ok {α β : Type} {x : KM α} {f : α → KM β} {b : β} (h : x >>= f = .ok b) :
    ∃ a, x = .ok a ∧ f a = .ok b := by
  cases x with
  | error e => cases h
  | ok a => exact ⟨a, rfl, h⟩

def guardK (c : Bool) (msg : String) : KM Unit := if c then .ok () else .error msg

lemma guardK_ok {c : Bool} {msg : String} {u : Unit} (h : guardK c msg = .ok u) : c = true := by
  unfold guardK at h; split at h
  · assumption
  · cases h

def checkBead (L : Link) (spec : BeadSpec) : KM Bead :=
  guardK (!(spec.ts.size = 0 || spec.ts.any (· ≥ L.n) || spec.start > 2 || spec.stop > 2))
    "malformed bead" >>= fun _ =>
  develop L spec.ts spec.steps >>= fun pos =>
  parallel (L.G[spec.ts[0]!]!) ((pos[0]!)[spec.start]!) (pos.back![spec.stop]!) >>= fun par =>
  guardK (!(par > 0)) "bead endpoints coincide" >>= fun _ =>
  (if par < 0 then piQ (L.G[spec.ts[0]!]!) ((pos[0]!)[spec.start]!) spec.piPoint
   else if spec.ts.size == 1 && spec.start == spec.stop then .error "degenerate one-triangle bead"
   else .ok none) >>= fun Q =>
  Bead.make spec.ts spec.start spec.stop (L.G[spec.ts[0]!]!) pos spec.steps Q >>= fun b =>
  checkCrossings (L.G[spec.ts[0]!]!) ((pos[0]!)[spec.start]!) b.outDir b.pi pos spec.steps
    (List.range (spec.ts.size - 1)) >>= fun _ =>
  .ok b

/-- An accepted bead: its geodesic, from `S` toward `outDir`, crosses every
interior edge of its developed strip (the minor arc to `E`, or for a π-bead
the half circle through `outDir`). -/
theorem checkBead_ok {L : Link} {spec : BeadSpec} {b : Bead} (h : checkBead L spec = .ok b) :
    b.S = (b.pos[0]!)[b.start]! ∧ b.E = b.pos.back![b.stop]! ∧
      (b.pi = false → b.outDir = b.E) ∧
      ∀ k ∈ List.range (b.ts.size - 1),
        crosses b.G b.S b.outDir ((b.pos[k]!)[(others3 b.steps[k]!).1]!)
          ((b.pos[k]!)[(others3 b.steps[k]!).2]!) b.pi = .ok true := by
  unfold checkBead at h
  obtain ⟨_, -, h⟩ := bind_ok h
  obtain ⟨pos, -, h⟩ := bind_ok h
  obtain ⟨par, -, h⟩ := bind_ok h
  obtain ⟨_, -, h⟩ := bind_ok h
  obtain ⟨Q, -, h⟩ := bind_ok h
  obtain ⟨b', hb, h⟩ := bind_ok h
  obtain ⟨_, hc, h⟩ := bind_ok h
  cases h
  obtain ⟨hts, hG, hposb, hst, hs, he, hS, hE, hout⟩ := Bead.make_ok hb
  refine ⟨by rw [hS, hposb, hs], by rw [hE, hposb, he], ?_, ?_⟩
  · intro hpi; rw [hout, hpi]; rfl
  · have := checkCrossings_ok hc
    rw [hts, hG, hposb, hst, hS]
    exact this

inductive Reason where
  | invalid (tet : Nat)
  | edge (e : Nat)
  | vertex (v : Nat) (necklace : Array BeadSpec)

/-- A certificate that the metric is not locally CAT(0) (`vlink-certificate/2`). -/
structure FalseCert where
  isosig : String
  sqlens : Array ℚ
  reason : Reason

def checkFalse (c : FalseCert) : KM String := do
  let T ← Tri.ofIsoSig c.isosig
  let geo ← Geometry.make T c.sqlens
  match c.reason with
  | .invalid t =>
    if t ≥ T.n || geo.tetValid t then throw s!"tetrahedron {t} is valid"
    return s!"tetrahedron {t} is degenerate"
  | r =>
    for t in [0:T.n] do
      if !geo.tetValid t then throw s!"tetrahedron {t} is degenerate (use reason \"invalid\")"
    match r with
    | .edge e =>
      if e ≥ T.edges.size then throw "bad edge"
      if atLeast (← geo.edgeAngle e) 2 then throw s!"edge {e} has dihedral sum >= 2*pi"
      return s!"edge {e} has dihedral sum < 2*pi"
    | .vertex iv necklace =>
      if iv ≥ T.corners.size then throw "bad vertex"
      let L ← Link.make geo iv
      let beads ← necklace.mapM (checkBead L)
      if beads.size = 0 then throw "empty necklace"
      if atLeast (← totalLength beads) 2 then throw "necklace length >= 2*pi"
      let m := beads.size
      for k in [0:m] do
        let b1 := beads[k]!
        let b2 := beads[(k + 1) % m]!
        if (L.lv[b1.ts.back!]!)[b1.stop]! ≠ (L.lv[b2.ts[0]!]!)[b2.start]! then
          throw s!"beads {k} and {(k + 1) % m} do not meet"
        for side in [1, 2] do
          if !atLeast (← junctionSide L b1 b2 side) 1 then
            throw s!"junction {k}: angle < pi on side {side}"
      return s!"vertex {iv}: necklace of {m} beads, length < 2*pi, all junctions >= pi"
    | .invalid _ => throw "unreachable"

/-! ### True certificates (docs/true_certificates.md) -/

def LEMMA_4_9 : String := "Elder-McCammond-Meier 2003, Lemma 4.9"

/-- For each pair `(x, y)` with `x` ending where `y` starts: does the junction
pass? Pairs involving a π-bead count as passing (conservative). -/
def junctionTable (L : Link) (beads : Array Bead) : KM (Std.HashMap (Nat × Nat) Bool) := do
  let mut byStart : Std.HashMap Nat (Array Nat) := {}
  for i in [0:beads.size] do
    let b := beads[i]!
    let w := (L.lv[b.ts[0]!]!)[b.start]!
    byStart := byStart.insert w ((byStart.getD w #[]).push i)
  let mut table := {}
  for x in [0:beads.size] do
    let bx := beads[x]!
    let w := (L.lv[bx.ts.back!]!)[bx.stop]!
    for y in byStart.getD w #[] do
      let by_ := beads[y]!
      let ok ← if bx.pi || by_.pi then pure true else junctionOk L bx by_
      table := table.insert (x, y) ok
  return table

/-- The chain `s .. x` of bead indices following predecessor pointers. -/
def chainOf (tree : Std.HashMap Nat Nat) (s x : Nat) : KM (Array Nat) := do
  let mut out := #[x]
  let mut x := x
  for _ in [0:tree.size + 2] do
    if x == s then return out.reverse
    match tree.get? x with
    | some p => x := p; out := out.push p
    | none => throw "chain does not reach the root"
  throw "predecessor pointers contain a cycle"

def checkVertexTrees (L : Link) (beads : Array Bead) (trees : Std.HashMap Nat (Std.HashMap Nat Nat)) :
    KM Unit := do
  let n := beads.size
  let lens := beads.map beadLength
  let J ← junctionTable L beads
  let mut succ : Std.HashMap Nat (Array Nat) := {}
  for ((x, y), ok) in J.toList do
    if ok then succ := succ.insert x ((succ.getD x #[]).push y)
  for s in [0:n] do
    let tree := trees.getD s {}
    for (x, _) in tree.toList do
      if x ≤ s || x ≥ n then throw s!"tree of {s}: bad bead index"
    -- exact chain sums, memoised per tree node
    let mut sums : Std.HashMap Nat AngleSum := {}
    sums := sums.insert s (← angleSum #[lens[s]!])
    let mut nodes := #[s]
    for (x, _) in tree.toList do
      let c ← chainOf tree s x
      for (u, w) in c.toList.zip c.toList.tail do
        if !(J.getD (u, w) false) then
          throw s!"tree of {s}: {u} -> {w} is not a passing junction"
      sums := sums.insert x (← angleSum (c.map (lens[·]!)))
      nodes := nodes.push x
    for x in nodes do
      let sx := sums.getD x AngleSum.zero
      for y in succ.getD x #[] do
        if y == s then
          if !atLeast sx 2 then
            throw s!"bead {s} closes a loop shorter than 2*pi through {x}"
          continue
        if y < s then continue
        let sxy ← angleAdd sx lens[y]!
        if atLeast sxy 2 then continue
        match sums.get? y with
        | none => throw s!"tree of {s}: {y} reachable via {x} but missing"
        | some sy =>
          if !(← angleLeq sy sxy) then
            throw s!"tree of {s}: distance to {y} not a lower bound via {x}"

/-- A certificate that the metric is locally CAT(0) (`vlink-certificate-true/2`). -/
structure TrueCert where
  isosig : String
  sqlens : Array ℚ
  piBeads : Bool
  allowRepeats : Bool
  assumptions : Array String
  /-- per vertex: number of beads and the trees -/
  vertices : Array (Nat × Nat × Std.HashMap Nat (Std.HashMap Nat Nat))

def checkTrue (c : TrueCert) : KM String := do
  if !c.piBeads && !c.assumptions.contains LEMMA_4_9 then
    throw s!"pi_beads=False requires the assumption: {LEMMA_4_9}"
  if !c.allowRepeats then throw "allow_repeats=False relies on an unproven rule"
  let T ← Tri.ofIsoSig c.isosig
  let geo ← Geometry.make T c.sqlens
  for t in [0:T.n] do
    if !geo.tetValid t then throw s!"tetrahedron {t} is degenerate"
  for e in [0:T.edges.size] do
    if !atLeast (← geo.edgeAngle e) 2 then throw s!"edge {e} has dihedral sum < 2*pi"
  let vs := c.vertices.qsort (fun a b => a.1 < b.1)
  if vs.map (·.1) != (List.range T.corners.size).toArray then
    throw "certificate must cover every vertex"
  for (iv, nBeads, trees) in vs do
    let L ← Link.make geo iv
    let beads ← enumerateBeads L c.piBeads
    if nBeads ≠ beads.size then
      throw s!"vertex {iv}: bead count differs from the kernel enumeration"
    checkVertexTrees L beads trees
  let note := if c.assumptions.isEmpty then "" else
    s!" (assuming {"; ".intercalate c.assumptions.toList})"
  return s!"locally CAT(0): {T.n} tetrahedra valid, all edges >= 2*pi, no short necklace at {vs.size} vertices{note}"

end VLink
