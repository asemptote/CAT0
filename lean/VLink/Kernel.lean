import VLink.Angle
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

abbrev KM := Except String

def kfuel : Nat := 64

def sgn (x : MQ) : KM Int :=
  match MQ.sign kfuel x with
  | some s => pure s
  | none => throw "sign undecided (out of fuel)"

def isZero (x : MQ) : KM Bool := return (← sgn x) == 0

def twSgn (x : Tw) : KM Int :=
  match Tw.sign kfuel x with
  | some s => pure s
  | none => throw "sign undecided (out of fuel)"

def mqOne : MQ := MQ.ofRat 1

/-! ### Vectors (unnormalised points of the link, coordinates in a fixed frame) -/

abbrev Vec := Array MQ
abbrev Gram := Array (Array ℚ)

def E3 : Array Vec := #[#[mqOne, [], []], #[[], mqOne, []], #[[], [], mqOne]]

def vsub (x y : Vec) : Vec := (x.zip y).map fun (a, b) => MQ.sub a b
def vneg (x : Vec) : Vec := x.map MQ.neg
def vscale (c : MQ) (x : Vec) : Vec := x.map (MQ.mul c)
def vscaleQ (c : ℚ) (x : Vec) : Vec := x.map (MQ.smul c)

def dot (G : Gram) (x y : Vec) : MQ := Id.run do
  let mut acc : MQ := []
  for j in [0:3] do
    for k in [0:3] do
      let g := (G[j]!)[k]!
      if g ≠ 0 then acc := MQ.add acc (MQ.smul g (MQ.mul x[j]! y[k]!))
  return acc

def det3 (x y z : Vec) : MQ :=
  let m := MQ.mul
  let s := MQ.sub
  MQ.add (s (m x[0]! (s (m y[1]! z[2]!) (m y[2]! z[1]!)))
            (m x[1]! (s (m y[0]! z[2]!) (m y[2]! z[0]!))))
         (m x[2]! (s (m y[0]! z[1]!) (m y[1]! z[0]!)))

/-- `+1` same direction, `-1` opposite, `0` not parallel. -/
def parallel (G : Gram) (x y : Vec) : KM Int := do
  let m := MQ.mul
  let c := #[MQ.sub (m x[1]! y[2]!) (m x[2]! y[1]!), MQ.sub (m x[2]! y[0]!) (m x[0]! y[2]!),
             MQ.sub (m x[0]! y[1]!) (m x[1]! y[0]!)]
  for v in c do
    if !(← isZero v) then return 0
  sgn (dot G x y)

/-! ### Angle sums -/

def angleAdd (st : AngleSum) (ab : MQ × MQ) : KM AngleSum :=
  match st.add ab.1 ab.2 with
  | some st' => pure st'
  | none => throw "angle undecided or between degenerate vectors"

def angleSum (l : Array (MQ × MQ)) : KM AngleSum :=
  l.foldlM angleAdd AngleSum.zero

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
  let normal := vsub (vsub Pd (vscaleQ tr.p.1 Pa)) (vscaleQ tr.p.2 Pb)
  let new := vsub (vscaleQ tr.q.1 Pa |> fun x => (x.zip (vscaleQ tr.q.2 Pb)).map
    fun (u, w) => MQ.add u w) (vscale tr.mu normal)
  let P2 := ((#[#[], #[], #[]] : Array Vec).set! tr.g[a]! Pa |>.set! tr.g[b]! Pb).set! tr.g[d]! new
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

def coneSigns (X A B C : Vec) : KM (Int × Int) := do
  let den ← sgn (det3 A B C)
  return ((← sgn (det3 X B C)) * den, (← sgn (det3 A X C)) * den)

/-- Does the geodesic from `S` toward `Q` cross the open edge arc `AB`? -/
def crosses (G : Gram) (S Q A B : Vec) (halfCircle : Bool) : KM Bool := do
  let dS := det3 A B S
  if ← isZero dS then return false
  let dQ := det3 A B Q
  let X := vsub (vscale dQ S) (vscale dS Q)
  let mut allZero := true
  for x in X do
    if !(← isZero x) then allZero := false
  if allZero then return false
  let (al, be) ← coneSigns X A B S
  if !(al == be && al != 0) then return false
  let Y := if al > 0 then X else vneg X
  if halfCircle then
    let qs := dot G Q S
    let ss := dot G S S
    return (← sgn (MQ.sub (MQ.mul (dot G Y Q) ss) (MQ.mul (dot G Y S) qs))) > 0
  let C := if !(← isZero (det3 S Q A)) then A else B
  let (ga, de) ← coneSigns Y S Q C
  return ga > 0 && de > 0

/-- `(a, b²)` of the angle at `W` between the geodesics toward `X` and `Y`. -/
def tangentAngle (G : Gram) (W X Y : Vec) : MQ × MQ :=
  let ww := dot G W W
  let xw := dot G X W
  let yw := dot G Y W
  let a := MQ.sub (MQ.mul (dot G X Y) ww) (MQ.mul xw yw)
  let nx := MQ.sub (MQ.mul (dot G X X) ww) (MQ.mul xw xw)
  let ny := MQ.sub (MQ.mul (dot G Y Y) ww) (MQ.mul yw yw)
  (a, MQ.sub (MQ.mul nx ny) (MQ.mul a a))

/-- `(a, b²)` of the angle between `S` and `E`. -/
def length (G : Gram) (S E : Vec) : MQ × MQ :=
  let a := dot G S E
  (a, MQ.sub (MQ.mul (dot G S S) (dot G E E)) (MQ.mul a a))

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
    (steps : Array Nat) (Q : Option Vec) : KM Bead := do
  let S := (pos[0]!)[start]!
  let E := pos.back![stop]!
  let pi := (← parallel G S E) < 0
  let q := Q.getD #[]
  return ⟨ts, start, stop, G, pos, steps, S, E, pi, if pi then q else S, if pi then q else E⟩

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
    let X : Vec := #[0, 1, 2].map fun i =>
      (List.range 3).foldl (fun acc j => MQ.add acc (MQ.mul bout.outDir[j]! (P[j]!)[i]!)) []
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

def checkBead (L : Link) (spec : BeadSpec) : KM Bead := do
  let ts := spec.ts
  let s := spec.start
  let e := spec.stop
  if ts.size = 0 || ts.any (· ≥ L.n) || s > 2 || e > 2 then throw "malformed bead"
  let pos ← develop L ts spec.steps
  let G := L.G[ts[0]!]!
  let S := (pos[0]!)[s]!
  let E := pos.back![e]!
  let par ← parallel G S E
  if par > 0 then throw "bead endpoints coincide"
  let mut Q : Option Vec := none
  if par < 0 then
    match spec.piPoint with
    | none => throw "pi-bead without pi_point"
    | some q =>
      if q.size ≠ 3 then throw "malformed pi_point"
      if (← parallel G S q) != 0 then throw "pi_point on the axis"
      Q := some q
  else if ts.size == 1 && s == e then throw "degenerate one-triangle bead"
  let b ← Bead.make ts s e G pos spec.steps Q
  for k in [0:ts.size - 1] do
    let (a, c) := others3 spec.steps[k]!
    let ok ← crosses G S (if b.pi then b.outDir else E) (pos[k]!)[a]! (pos[k]!)[c]! b.pi
    if !ok then throw s!"geodesic of bead {ts} does not cross interior edge {k}"
  return b

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
