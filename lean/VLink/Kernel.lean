import VLink.Develop
import VLink.Tri
import VLink.Necklace

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

/-! ### Lists in `KM` -/

/-- `f` succeeds on every element of the list. -/
def forallK {α : Type} : List α → (α → KM Unit) → KM Unit
  | [], _ => .ok ()
  | a :: l, f =>
    match f a with
    | .ok () => forallK l f
    | .error e => .error e

lemma forallK_ok {α : Type} {l : List α} {f : α → KM Unit} (h : forallK l f = .ok ()) :
    ∀ a ∈ l, f a = .ok () := by
  induction l with
  | nil => simp
  | cons a l ih =>
    unfold forallK at h
    split at h
    · rename_i ha
      intro b hb
      rcases List.mem_cons.mp hb with rfl | hb
      · exact ha
      · exact ih h b hb
    · cases h

/-- `List.mapM` in `KM`, with a usable specification. -/
def mapK {α β : Type} (f : α → KM β) : List α → KM (List β)
  | [] => .ok []
  | a :: l =>
    match f a, mapK f l with
    | .ok b, .ok bs => .ok (b :: bs)
    | .error e, _ => .error e
    | _, .error e => .error e

lemma mapK_ok {α β : Type} {f : α → KM β} {l : List α} {l' : List β}
    (h : mapK f l = .ok l') : List.Forall₂ (fun a b => f a = .ok b) l l' := by
  induction l generalizing l' with
  | nil => cases h; exact List.Forall₂.nil
  | cons a l ih =>
    unfold mapK at h
    split at h
    · rename_i b bs ha hl
      cases h
      exact List.Forall₂.cons ha (ih hl)
    · cases h
    · cases h

lemma forall2_mem_right {α β : Type} {R : α → β → Prop} {l : List α} {l' : List β}
    (h : List.Forall₂ R l l') {b : β} (hb : b ∈ l') : ∃ a ∈ l, R a b := by
  induction h with
  | nil => simp at hb
  | cons hab _ ih =>
    rcases List.mem_cons.mp hb with rfl | hb
    · exact ⟨_, List.mem_cons_self, hab⟩
    · obtain ⟨a, ha, hR⟩ := ih hb
      exact ⟨a, List.mem_cons_of_mem _ ha, hR⟩

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

/-- `(a, b²)` of the dihedral angles around edge `e`. -/
def Geometry.edgeAngles (g : Geometry) (e : Nat) : Array (MQ × MQ) :=
  g.T.edges[e]!.map fun (tet, k) => g.dihedral tet EDGES[k]!.1 EDGES[k]!.2

def Geometry.edgeAngle (g : Geometry) (e : Nat) : KM AngleSum := angleSum (g.edgeAngles e)

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
  trans : List (List Trans)
  deriving Inhabited

def others3 (d : Nat) : Nat × Nat := if d = 0 then (1, 2) else if d = 1 then (0, 2) else (0, 1)

lemma others3_spec {d : Nat} (hd : d < 3) :
    (others3 d).1 < 3 ∧ (others3 d).2 < 3 ∧ (others3 d).1 ≠ (others3 d).2 ∧
      (others3 d).1 ≠ d ∧ (others3 d).2 ≠ d := by
  unfold others3; split_ifs <;> omega

/-- Entry `(i, j)` of a Gram matrix. -/
def gq (G : Gram) (i j : Nat) : ℚ := (G[i]!)[j]!

/-- `g` permutes `{0, 1, 2}`. -/
def isPerm3 (g : Array Nat) : Bool :=
  g[0]! < 3 && g[1]! < 3 && g[2]! < 3 && g[0]! != g[1]! && g[0]! != g[2]! && g[1]! != g[2]!

def symm3 (G : Gram) : Bool := gq G 0 1 == gq G 1 0 && gq G 0 2 == gq G 2 0 && gq G 1 2 == gq G 2 1

lemma isPerm3_spec {g : Array Nat} (h : isPerm3 g = true) {i j : Nat} (hi : i < 3) (hj : j < 3)
    (hij : i ≠ j) : g[i]! < 3 ∧ g[i]! ≠ g[j]! := by
  unfold isPerm3 at h
  simp only [Bool.and_eq_true, decide_eq_true_eq, bne_iff_ne, ne_eq] at h
  interval_cases i <;> interval_cases j <;> omega

lemma symm3_spec {G : Gram} (h : symm3 G = true) {i j : Nat} (hi : i < 3) (hj : j < 3) :
    gq G i j = gq G j i := by
  unfold symm3 at h
  simp only [Bool.and_eq_true, beq_iff_eq] at h
  interval_cases i <;> interval_cases j <;> simp_all

lemma bform_symm {G : Gram} (h : symm3 G = true) (u v : Fin 3 → ℝ) :
    bform G.toMat u v = bform G.toMat v u := by
  have s := fun i j (hi : i < 3) (hj : j < 3) => symm3_spec h hi hj
  simp only [bform, Fin.sum_univ_three, Gram.toMat]
  have e01 := s 0 1 (by norm_num) (by norm_num)
  have e02 := s 0 2 (by norm_num) (by norm_num)
  have e12 := s 1 2 (by norm_num) (by norm_num)
  simp only [gq] at e01 e02 e12
  simp only [Fin.val_zero, Fin.val_one, Fin.val_two, e01, e02, e12]
  ring

/-- The leading principal minors of `G` are positive (Sylvester's criterion). -/
def posdef3 (G : Gram) : Bool :=
  0 < gq G 0 0 && 0 < gq G 0 0 * gq G 1 1 - gq G 0 1 ^ 2 &&
    0 < gq G 0 0 * gq G 1 1 * gq G 2 2 + 2 * gq G 0 1 * gq G 0 2 * gq G 1 2 -
      gq G 0 0 * gq G 1 2 ^ 2 - gq G 0 1 ^ 2 * gq G 2 2 - gq G 0 2 ^ 2 * gq G 1 1

/-- Sylvester's criterion for a symmetric `3 × 3` matrix: positive leading
minors make the form positive definite. -/
theorem posdef3_spec {G : Gram} (hs : symm3 G = true) (hp : posdef3 G = true)
    (x : Fin 3 → ℝ) (hx : x ≠ 0) : 0 < bform G.toMat x x := by
  unfold posdef3 at hp
  simp only [Bool.and_eq_true, decide_eq_true_eq] at hp
  obtain ⟨⟨ha, hm⟩, hD⟩ := hp
  have s01 := symm3_spec hs (i := 1) (j := 0) (by norm_num) (by norm_num)
  have s02 := symm3_spec hs (i := 2) (j := 0) (by norm_num) (by norm_num)
  have s12 := symm3_spec hs (i := 2) (j := 1) (by norm_num) (by norm_num)
  set a : ℝ := ((gq G 0 0 : ℚ) : ℝ)
  set b : ℝ := ((gq G 0 1 : ℚ) : ℝ)
  set c : ℝ := ((gq G 0 2 : ℚ) : ℝ)
  set d : ℝ := ((gq G 1 1 : ℚ) : ℝ)
  set e : ℝ := ((gq G 1 2 : ℚ) : ℝ)
  set f : ℝ := ((gq G 2 2 : ℚ) : ℝ)
  have ha' : 0 < a := Rat.cast_pos.mpr ha
  have hm' : 0 < a * d - b ^ 2 := by
    have := (Rat.cast_pos (K := ℝ)).mpr hm; push_cast at this; exact this
  have hD' : 0 < a * d * f + 2 * b * c * e - a * e ^ 2 - b ^ 2 * f - c ^ 2 * d := by
    have := (Rat.cast_pos (K := ℝ)).mpr hD; push_cast at this; exact this
  have hQ : bform G.toMat x x = a * x 0 ^ 2 + d * x 1 ^ 2 + f * x 2 ^ 2 + 2 * b * x 0 * x 1 +
      2 * c * x 0 * x 2 + 2 * e * x 1 * x 2 := by
    have e10 : (((G[1]!)[0]! : ℚ) : ℝ) = b := by simp only [b, gq] at s01 ⊢; rw [s01]
    have e20 : (((G[2]!)[0]! : ℚ) : ℝ) = c := by simp only [c, gq] at s02 ⊢; rw [s02]
    have e21 : (((G[2]!)[1]! : ℚ) : ℝ) = e := by simp only [e, gq] at s12 ⊢; rw [s12]
    simp only [bform, Fin.sum_univ_three, Gram.toMat, Fin.val_zero, Fin.val_one, Fin.val_two, e10,
      e20, e21]
    simp only [a, b, c, d, e, f, gq]
    ring
  set m := a * d - b ^ 2
  set D := a * d * f + 2 * b * c * e - a * e ^ 2 - b ^ 2 * f - c ^ 2 * d
  have key : a * m * bform G.toMat x x = m * (a * x 0 + b * x 1 + c * x 2) ^ 2 +
      (m * x 1 + (a * e - b * c) * x 2) ^ 2 + a * D * x 2 ^ 2 := by
    rw [hQ]; simp only [m, D]; ring
  have hR : 0 < m * (a * x 0 + b * x 1 + c * x 2) ^ 2 + (m * x 1 + (a * e - b * c) * x 2) ^ 2 +
      a * D * x 2 ^ 2 := by
    have t1 : 0 ≤ m * (a * x 0 + b * x 1 + c * x 2) ^ 2 := by positivity
    have t2 : 0 ≤ (m * x 1 + (a * e - b * c) * x 2) ^ 2 := by positivity
    by_cases h2 : x 2 = 0
    · by_cases h1 : x 1 = 0
      · have h0 : x 0 ≠ 0 := by
          intro h0; apply hx; ext i; fin_cases i <;> simp [h0, h1, h2]
        have : 0 < m * (a * x 0) ^ 2 := by
          have : 0 < (a * x 0) ^ 2 := by positivity
          positivity
        rw [h1, h2]; simp only [mul_zero, add_zero, zero_pow two_ne_zero, mul_zero] ; nlinarith
      · have : 0 < (m * x 1) ^ 2 := by positivity
        rw [h2]; simp only [mul_zero, add_zero, zero_pow two_ne_zero]; nlinarith
    · have : 0 < a * D * x 2 ^ 2 := by positivity
      linarith
  have ham : 0 < a * m := mul_pos ha' hm'
  rw [← key] at hR
  exact pos_of_mul_pos_right hR ham.le

def mkTrans (G : Array Gram) (adj : Array (Array Nat)) (glue : Array (Array (Array Nat)))
    (t d : Nat) : KM Trans :=
  let t2 := (adj[t]!)[d]!
  let g := (glue[t]!)[d]!
  let a := (others3 d).1
  let b := (others3 d).2
  let Gt := G[t]!
  let G2 := G[t2]!
  if !(isPerm3 g && symm3 Gt && symm3 G2 && posdef3 Gt && posdef3 G2) then .error "malformed link data"
  else if gq Gt a a ≠ gq G2 g[a]! g[a]! ∨ gq Gt a b ≠ gq G2 g[a]! g[b]! ∨
      gq Gt b b ≠ gq G2 g[b]! g[b]! then
    .error "glued link edge has different lengths on its two sides"
  else
    let dt := gq Gt a a * gq Gt b b - gq Gt a b ^ 2
    let p0 := (gq Gt d a * gq Gt b b - gq Gt d b * gq Gt a b) / dt
    let p1 := (gq Gt d b * gq Gt a a - gq Gt d a * gq Gt a b) / dt
    let q0 := (gq G2 g[d]! g[a]! * gq Gt b b - gq G2 g[d]! g[b]! * gq Gt a b) / dt
    let q1 := (gq G2 g[d]! g[b]! * gq Gt a a - gq G2 g[d]! g[a]! * gq Gt a b) / dt
    let h := gq Gt d d - p0 * gq Gt d a - p1 * gq Gt d b
    let h2 := gq G2 g[d]! g[d]! - q0 * gq G2 g[d]! g[a]! - q1 * gq G2 g[d]! g[b]!
    if dt = 0 ∨ h ≤ 0 ∨ h2 ≤ 0 then .error "degenerate link triangle"
    else .ok ⟨t2, g, (p0, p1), (q0, q1), MQ.sqrtRat (h2 / h)⟩

/-- The link's gluing is consistent: crossing from `t` over the edge opposite
`d` lands in a link triangle and crosses back to `t` with the inverse vertex
map, and the two vertices of the shared edge keep their link-vertex labels. -/
def glueOk (L : Link) : Bool :=
  (List.range L.n).all fun t => (List.range 3).all fun d =>
    let t2 := (L.adj[t]!)[d]!
    let g := (L.glue[t]!)[d]!
    t2 < L.n && isPerm3 g && (L.adj[t2]!)[g[d]!]! == t &&
      (List.range 3).all (fun j => ((L.glue[t2]!)[g[d]!]!)[g[j]!]! == j) &&
      (List.range 3).all (fun j => j == d || (L.lv[t2]!)[g[j]!]! == (L.lv[t]!)[j]!)

lemma glueOk_spec {L : Link} (h : glueOk L = true) {t d : Nat} (ht : t < L.n) (hd : d < 3) :
    (L.adj[t]!)[d]! < L.n ∧ isPerm3 ((L.glue[t]!)[d]!) = true ∧
      (L.adj[(L.adj[t]!)[d]!]!)[((L.glue[t]!)[d]!)[d]!]! = t ∧
      (∀ j < 3, j ≠ d → (L.lv[(L.adj[t]!)[d]!]!)[((L.glue[t]!)[d]!)[j]!]! = (L.lv[t]!)[j]!) := by
  unfold glueOk at h
  simp only [List.all_eq_true, List.mem_range, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq,
    Bool.or_eq_true] at h
  obtain ⟨⟨⟨⟨h1, h2⟩, h3⟩, -⟩, h5⟩ := h t ht d hd
  exact ⟨h1, h2, h3, fun j hj hjd => (h5 j hj).resolve_left hjd⟩

/-- `c, a, b` are the three vertices of a link triangle, in some order. -/
def Corner3 (c a b : Nat) : Prop := c < 3 ∧ a < 3 ∧ b < 3 ∧ c ≠ a ∧ c ≠ b ∧ a ≠ b

/-- One step of the fan walk around a link vertex: from triangle `t`, whose
vertex `c` is at the link vertex, cross the edge opposite `ent` (the edge
from `c` to `ex`) into the adjacent triangle. -/
def fanNext (L : Link) (s : Nat × Nat × Nat × Nat) : Nat × Nat × Nat × Nat :=
  let g := (L.glue[s.1]!)[s.2.2.1]!
  ((L.adj[s.1]!)[s.2.2.1]!, g[s.2.1]!, g[s.2.2.2]!, 3 - g[s.2.1]! - g[s.2.2.2]!)

/-- One step of the fan walk stays at the same link vertex and moves into the
triangle glued to `t` along the edge from it to the vertex `ex`; in the new
triangle that edge runs from the link vertex to the new `ent`, and the new
triangle is glued back to `t` along it. So the walk turns around the link
vertex one triangle at a time. -/
theorem fan_step {L : Link} (hg : glueOk L = true) {t c ent ex : Nat} (ht : t < L.n)
    (hc : Corner3 c ent ex) :
    (fanNext L (t, c, ent, ex)).1 < L.n ∧
      Corner3 (fanNext L (t, c, ent, ex)).2.1 (fanNext L (t, c, ent, ex)).2.2.1
        (fanNext L (t, c, ent, ex)).2.2.2 ∧
      (L.lv[(fanNext L (t, c, ent, ex)).1]!)[(fanNext L (t, c, ent, ex)).2.1]! = (L.lv[t]!)[c]! ∧
      (L.lv[(fanNext L (t, c, ent, ex)).1]!)[(fanNext L (t, c, ent, ex)).2.2.1]! = (L.lv[t]!)[ex]! ∧
      (L.adj[(fanNext L (t, c, ent, ex)).1]!)[((L.glue[t]!)[ent]!)[ent]!]! = t := by
  obtain ⟨hc3, he3, hx3, hce, hcx, hex⟩ := hc
  obtain ⟨ht2, hperm, hback, hlv⟩ := glueOk_spec hg ht he3
  obtain ⟨gc3, gce⟩ := isPerm3_spec hperm hc3 he3 hce
  obtain ⟨-, gcx⟩ := isPerm3_spec hperm hc3 hx3 hcx
  obtain ⟨gx3, gxe⟩ := isPerm3_spec hperm hx3 he3 (Ne.symm hex)
  simp only [fanNext]
  exact ⟨ht2, ⟨gc3, gx3, by omega, gcx, by omega, by omega⟩, hlv c hc3 hce, hlv ex hx3 (Ne.symm hex), hback⟩


/-- Link data without the transitions. -/
def Link.base (geo : Geometry) (iv : Nat) : Link := Id.run do
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
  return ⟨corners.size, G, adj, glue, lv, labels.size, []⟩

def Link.make (geo : Geometry) (iv : Nat) : KM Link :=
  let B := Link.base geo iv
  if !glueOk B then .error "inconsistent link gluing" else
  match mapK (fun t => mapK (mkTrans B.G B.adj B.glue t) [0, 1, 2]) (List.range B.n) with
  | .ok rows => .ok { B with trans := rows }
  | .error e => .error e

theorem Link.make_glueOk {geo : Geometry} {iv : Nat} {L : Link} (h : Link.make geo iv = .ok L) :
    glueOk L = true := by
  unfold Link.make at h
  simp only at h
  split_ifs at h with hgl
  split at h
  · cases h
    simp only [Bool.not_eq_true', Bool.not_eq_false] at hgl
    exact hgl
  · cases h

/-- The transition across the edge opposite vertex `d` of triangle `t`. -/
def Link.tr (L : Link) (t d : Nat) : Trans := (L.trans.getD t []).getD d default

theorem Link.make_trans {geo : Geometry} {iv : Nat} {L : Link} (h : Link.make geo iv = .ok L) :
    ∀ t < L.n, ∀ d < 3, mkTrans L.G L.adj L.glue t d = .ok (L.tr t d) := by
  unfold Link.make at h
  simp only at h
  split_ifs at h with hgl
  split at h
  · rename_i rows hrows
    cases h
    intro t ht d hd
    dsimp only at ht
    have hf := mapK_ok hrows
    have hlen : rows.length = (Link.base geo iv).n := by simpa using hf.length_eq.symm
    have htr : t < rows.length := by omega
    have ht' := (List.forall₂_iff_get.mp hf).2 t (by simpa using ht) htr
    simp only [List.get_eq_getElem, List.getElem_range] at ht'
    have hf2 := mapK_ok ht'
    have hlen2 : rows[t].length = 3 := by simpa using hf2.length_eq.symm
    have hd' := (List.forall₂_iff_get.mp hf2).2 d (by simpa using hd) (by omega)
    simp only [Link.tr]
    rw [List.getD_eq_getElem _ _ htr, List.getD_eq_getElem _ _ (by omega)]
    interval_cases d <;> simpa using hd'
  · cases h

/-- Positions of a link triangle's vertices, developed into the frame of
Gram matrix `M`: `stepWith` places the neighbour's vertices. -/
def stepWith (tr : Trans) (P : Array V3) (d : Nat) : Array V3 :=
  let a := (others3 d).1
  let b := (others3 d).2
  let Pa := P[a]!
  let Pb := P[b]!
  let Pd := P[d]!
  let normal := V3.sub (V3.sub Pd (V3.scaleQ tr.p.1 Pa)) (V3.scaleQ tr.p.2 Pb)
  let new := V3.sub (V3.add (V3.scaleQ tr.q.1 Pa) (V3.scaleQ tr.q.2 Pb)) (V3.scale tr.mu normal)
  let f : Nat → V3 := fun i => if i = tr.g[a]! then Pa else if i = tr.g[b]! then Pb else new
  #[f 0, f 1, f 2]

/-- Develop across the edge opposite vertex `d` of triangle `t` at positions `P`. -/
def step (L : Link) (P : Array Vec) (t d : Nat) : Nat × Array Vec :=
  ((L.tr t d).t2, stepWith (L.tr t d) P d)

def developFrom (L : Link) : Array V3 → List Nat → List Nat → KM (List (Array V3))
  | _, [_], [] => .ok []
  | P, t :: t' :: ts, d :: ds =>
    if d > 2 ∨ (L.adj[t]!)[d]! ≠ t' then .error s!"strip step does not glue triangle {t} to {t'}"
    else
      match developFrom L (step L P t d).2 (t' :: ts) ds with
      | .ok rest => .ok ((step L P t d).2 :: rest)
      | .error e => .error e
  | _, _, _ => .error "wrong number of strip steps"

def develop (L : Link) (ts steps : Array Nat) : KM (Array (Array Vec)) :=
  if steps.size + 1 ≠ ts.size then .error "wrong number of strip steps"
  else match developFrom L E3 ts.toList steps.toList with
    | .ok rest => .ok (E3 :: rest).toArray
    | .error e => .error e

/-- The vertices `P` of a triangle with Gram matrix `Gt`, in the frame `M`. -/
def GramOK (M : Fin 3 → Fin 3 → ℝ) (P : Array V3) (Gt : Gram) : Prop :=
  ∀ i j, i < 3 → j < 3 → bform M (P[i]!).val (P[j]!).val = (gq Gt i j : ℝ)

/-- One development step is isometric: the neighbour's developed vertices
have its Gram matrix, and its new vertex lies on the other side of the
shared edge's plane. -/
theorem stepWith_iso {G : Array Gram} {adj : Array (Array Nat)} {glue : Array (Array (Array Nat))}
    {t d : Nat} {tr : Trans} {M : Fin 3 → Fin 3 → ℝ} (hsym : ∀ u v, bform M u v = bform M v u)
    (h : mkTrans G adj glue t d = .ok tr) (hd : d < 3) {P : Array V3} (hP : GramOK M P G[t]!) :
    GramOK M (stepWith tr P d) G[tr.t2]! ∧
      ∃ μ : ℝ, 0 < μ ∧ rdet (P[(others3 d).1]!).val (P[(others3 d).2]!).val
        ((stepWith tr P d)[tr.g[d]!]!).val =
        -μ * rdet (P[(others3 d).1]!).val (P[(others3 d).2]!).val (P[d]!).val := by
  obtain ⟨ha3, hb3, hab, had, hbd⟩ := others3_spec hd
  unfold mkTrans at h
  simp only at h
  split_ifs at h with h1 h2 h3 <;> try (cases h; done)
  cases h
  set a := (others3 d).1 with ha
  set b := (others3 d).2 with hb
  set g := (glue[t]!)[d]! with hg
  set t2 := (adj[t]!)[d]! with ht2
  set Gt := G[t]! with hGt
  set G2 := G[t2]! with hG2
  simp only [Bool.not_eq_true', Bool.and_eq_false_iff, Bool.not_eq_false, not_or] at h1
  have hperm : isPerm3 g = true := h1.1.1.1.1
  have hsGt : symm3 Gt = true := h1.1.1.1.2
  have hsG2 : symm3 G2 = true := h1.1.1.2
  push Not at h2 h3
  obtain ⟨eaa, eab, ebb⟩ := h2
  obtain ⟨hdt, hh, hh2⟩ := h3
  -- the permutation g
  obtain ⟨gab3, gab⟩ := isPerm3_spec hperm ha3 hb3 hab
  obtain ⟨-, gad⟩ := isPerm3_spec hperm ha3 hd had
  obtain ⟨gb3, gbd⟩ := isPerm3_spec hperm hb3 hd hbd
  obtain ⟨gd3, -⟩ := isPerm3_spec hperm hd ha3 (Ne.symm had)
  -- real data
  set Pa := (P[a]!).val
  set Pb := (P[b]!).val
  set Pd := (P[d]!).val
  set p0 := (gq Gt d a * gq Gt b b - gq Gt d b * gq Gt a b) / (gq Gt a a * gq Gt b b - gq Gt a b ^ 2)
  set p1 := (gq Gt d b * gq Gt a a - gq Gt d a * gq Gt a b) / (gq Gt a a * gq Gt b b - gq Gt a b ^ 2)
  set q0 := (gq G2 g[d]! g[a]! * gq Gt b b - gq G2 g[d]! g[b]! * gq Gt a b) /
    (gq Gt a a * gq Gt b b - gq Gt a b ^ 2)
  set q1 := (gq G2 g[d]! g[b]! * gq Gt a a - gq G2 g[d]! g[a]! * gq Gt a b) /
    (gq Gt a a * gq Gt b b - gq Gt a b ^ 2)
  set hq := gq Gt d d - p0 * gq Gt d a - p1 * gq Gt d b
  set h2q := gq G2 g[d]! g[d]! - q0 * gq G2 g[d]! g[a]! - q1 * gq G2 g[d]! g[b]!
  set μ := (MQ.sqrtRat (h2q / hq)).val with hμdef
  have hμ0 : 0 < μ := by
    rw [hμdef, MQ.val_sqrtRat _ (le_of_lt (div_pos hh2 hh))]
    exact Real.sqrt_pos.mpr (by exact_mod_cast div_pos hh2 hh)
  have hμ2 : μ ^ 2 * (hq : ℝ) = (h2q : ℝ) := by
    rw [hμdef, MQ.val_sqrtRat _ (le_of_lt (div_pos hh2 hh)),
      Real.sq_sqrt (by exact_mod_cast (div_pos hh2 hh).le)]
    have : (hq : ℝ) ≠ 0 := by exact_mod_cast hh.ne'
    push_cast
    field_simp
  obtain ⟨hna, hnb, hnn, hdet⟩ := step_iso M hsym Pa Pb Pd
    (gq Gt a a) (gq Gt a b) (gq Gt b b) (gq Gt d a) (gq Gt d b) (gq Gt d d)
    (gq G2 g[d]! g[a]!) (gq G2 g[d]! g[b]!) (gq G2 g[d]! g[d]!) p0 p1 q0 q1 μ
    (hP a a ha3 ha3) (hP a b ha3 hb3) (hP b b hb3 hb3) (hP d a hd ha3) (hP d b hd hb3)
    (hP d d hd hd) (by exact_mod_cast hdt) (by simp [p0]) (by simp [p1]) (by simp [q0])
    (by simp [q1]) (by have e := hμ2; simp only [hq, h2q, Rat.cast_sub, Rat.cast_mul] at e; exact e)
  -- the developed positions
  set newR := (q0 : ℝ) • Pa + (q1 : ℝ) • Pb - μ • (Pd - (p0 : ℝ) • Pa - (p1 : ℝ) • Pb)
  have hE : ∀ i < 3, (stepWith ⟨t2, g, (p0, p1), (q0, q1), MQ.sqrtRat (h2q / hq)⟩ P d)[i]! =
      (if i = g[a]! then P[a]! else if i = g[b]! then P[b]! else
        V3.sub (V3.add (V3.scaleQ q0 P[a]!) (V3.scaleQ q1 P[b]!))
          (V3.scale (MQ.sqrtRat (h2q / hq))
            (V3.sub (V3.sub P[d]! (V3.scaleQ p0 P[a]!)) (V3.scaleQ p1 P[b]!)))) := by
    intro i hi
    interval_cases i <;> rfl
  have vA : (stepWith ⟨t2, g, (p0, p1), (q0, q1), MQ.sqrtRat (h2q / hq)⟩ P d)[g[a]!]!.val = Pa := by
    rw [hE _ gab3, if_pos rfl]
  have vB : (stepWith ⟨t2, g, (p0, p1), (q0, q1), MQ.sqrtRat (h2q / hq)⟩ P d)[g[b]!]!.val = Pb := by
    rw [hE _ gb3, if_neg (Ne.symm gab), if_pos rfl]
  have vD : (stepWith ⟨t2, g, (p0, p1), (q0, q1), MQ.sqrtRat (h2q / hq)⟩ P d)[g[d]!]!.val = newR := by
    rw [hE _ gd3, if_neg (Ne.symm gad), if_neg (Ne.symm gbd)]
    simp [newR, Pa, Pb, Pd, μ]
  have sG2 := fun i j (hi : i < 3) (hj : j < 3) => symm3_spec hsG2 hi hj
  have sGt := fun i j (hi : i < 3) (hj : j < 3) => symm3_spec hsGt hi hj
  refine ⟨?_, μ, hμ0, ?_⟩
  · intro i j hi hj
    have cover : ∀ k < 3, k = g[a]! ∨ k = g[b]! ∨ k = g[d]! := by intro k hk; omega
    rcases cover i hi with rfl | rfl | rfl <;> rcases cover j hj with rfl | rfl | rfl
    · rw [vA, hP a a ha3 ha3, eaa]
    · rw [vA, vB, hP a b ha3 hb3, eab]
    · rw [vA, vD, hsym, hna, sG2 _ _ gd3 gab3]
    · rw [vB, vA, hP b a hb3 ha3, sGt _ _ hb3 ha3, eab, sG2 _ _ gab3 gb3]
    · rw [vB, hP b b hb3 hb3, ebb]
    · rw [vB, vD, hsym, hnb, sG2 _ _ gd3 gb3]
    · rw [vD, vA, hna]
    · rw [vD, vB, hnb]
    · rw [vD, hnn]
  · rw [vD]
    exact hdet



lemma mkTrans_t2 {G : Array Gram} {adj : Array (Array Nat)} {glue : Array (Array (Array Nat))}
    {t d : Nat} {tr : Trans} (h : mkTrans G adj glue t d = .ok tr) :
    tr.t2 = (adj[t]!)[d]! ∧ symm3 G[t]! = true ∧ posdef3 G[t]! = true := by
  unfold mkTrans at h
  simp only at h
  split_ifs at h with h1 h2 h3 <;> try (cases h; done)
  cases h
  simp only [Bool.not_eq_true', Bool.and_eq_false_iff, Bool.not_eq_false, not_or] at h1
  exact ⟨rfl, h1.1.1.1.2, h1.1.2⟩

/-- Every link triangle's Gram matrix is positive definite: the developed
vectors live in a Euclidean space. -/
theorem link_gram_posdef {geo : Geometry} {iv : Nat} {L : Link} (hL : Link.make geo iv = .ok L)
    {t : Nat} (ht : t < L.n) (x : Fin 3 → ℝ) (hx : x ≠ 0) : 0 < bform (L.G[t]!).toMat x x := by
  obtain ⟨-, hs, hp⟩ := mkTrans_t2 (Link.make_trans hL t ht 0 (by norm_num))
  exact posdef3_spec hs hp x hx

lemma E3_gram (G0 : Gram) : GramOK G0.toMat E3 G0 := by
  intro i j hi hj
  interval_cases i <;> interval_cases j <;>
    simp [bform, E3, V3.val, Fin.sum_univ_three, Gram.toMat, gq, mqOne]

/-- Developing along a strip keeps every triangle's Gram matrix. -/
theorem developFrom_iso {L : Link}
    (hL : ∀ t < L.n, ∀ d < 3, mkTrans L.G L.adj L.glue t d = .ok (L.tr t d))
    {M : Fin 3 → Fin 3 → ℝ} (hsym : ∀ u v, bform M u v = bform M v u) :
    ∀ (ds ts : List Nat) (t : Nat) (P : Array V3) (out : List (Array V3)),
      t < L.n → (∀ x ∈ ts, x < L.n) → GramOK M P L.G[t]! →
      developFrom L P (t :: ts) ds = .ok out →
      ∀ k (hk : k < out.length), GramOK M out[k] L.G[ts[k]!]! := by
  intro ds
  induction ds with
  | nil =>
    intro ts t P out _ _ _ h k hk
    cases ts with
    | nil => simp [developFrom] at h; subst h; simp at hk
    | cons t' ts => simp [developFrom] at h
  | cons d ds ih =>
    intro ts t P out ht hts hP h k hk
    cases ts with
    | nil => simp [developFrom] at h
    | cons t' ts =>
      unfold developFrom at h
      split_ifs at h with hbad
      push Not at hbad
      obtain ⟨hd, hadj⟩ := hbad
      have hd3 : d < 3 := by omega
      have hm := hL t ht d hd3
      obtain ⟨ht2, -, -⟩ := mkTrans_t2 hm
      have hstep := (stepWith_iso hsym hm hd3 hP).1
      rw [ht2, hadj] at hstep
      split at h
      · rename_i rest hrest
        cases h
        have ht' : t' < L.n := hts t' List.mem_cons_self
        have hts' : ∀ x ∈ ts, x < L.n := fun x hx => hts x (List.mem_cons_of_mem _ hx)
        cases k with
        | zero => simpa [step] using hstep
        | succ k =>
          have := ih ts t' (step L P t d).2 rest ht' hts' (by simpa [step] using hstep) hrest k
            (by simpa using hk)
          simpa using this
      · cases h

/-- A developed strip is isometric: in the frame of its first triangle, the
`k`-th developed triangle's vertex vectors have exactly the Gram matrix of
the strip's `k`-th link triangle. -/
theorem develop_iso {geo : Geometry} {iv : Nat} {L : Link} (hL : Link.make geo iv = .ok L)
    {ts steps : Array Nat} {pos : Array (Array V3)} (h : develop L ts steps = .ok pos)
    (hts : ∀ t ∈ ts.toList, t < L.n) :
    ∀ k (hk : k < pos.size), GramOK (L.G[ts[0]!]!).toMat pos[k] L.G[ts[k]!]! := by
  have hLt := Link.make_trans hL
  unfold develop at h
  split_ifs at h with hsz
  split at h
  · rename_i rest hrest
    cases h
    obtain ⟨t0, tl, htl⟩ : ∃ t0 tl, ts.toList = t0 :: tl := by
      cases hc : ts.toList with
      | nil =>
        have : ts.size = 0 := by simpa using congrArg List.length hc
        omega
      | cons t0 tl => exact ⟨t0, tl, rfl⟩
    have h0 : ts[0]! = t0 := by
      have h' : ts.toList[0]! = t0 := by rw [htl]; rfl
      rwa [Array.getElem!_toList] at h'
    have ht0 : t0 < L.n := hts t0 (by rw [htl]; exact List.mem_cons_self)
    have hsym := bform_symm (mkTrans_t2 (hLt t0 ht0 0 (by norm_num))).2.1
    rw [htl] at hrest
    intro k hk
    rw [h0]
    cases k with
    | zero => simpa [h0] using E3_gram L.G[t0]!
    | succ k =>
      have := developFrom_iso hLt hsym steps.toList tl t0 E3 rest ht0
        (fun x hx => hts x (by rw [htl]; exact List.mem_cons_of_mem _ hx)) (E3_gram _) hrest k
        (by simpa using hk)
      have hk' : ts[k + 1]! = tl[k]! := by
        have h' : ts.toList[k + 1]! = tl[k]! := by rw [htl]; rfl
        rwa [Array.getElem!_toList] at h'
      simpa [hk'] using this
  · cases h

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

/-! ### Bead enumeration (completeness: `VLink/Complete.lean`, docs/proofs.md §3) -/

/-- The directions at `S` toward the open arc `AB`, as an oriented pair
`(L, R)` with `det(S, L, R) > 0`; `none` if `S` is on the great circle of `AB`. -/
def visible (S A B : Vec) : KM (Option (Vec × Vec)) :=
  match sgn (det3 S A B) with
  | .error e => .error e
  | .ok o => .ok (if o = 0 then none else some (if o > 0 then (A, B) else (B, A)))

/-- Is `X` inside the interval of directions `(L, R)` at `S`? -/
def inside (S L R X : Vec) (strict : Bool) : KM Bool :=
  match sgn (det3 S L X), sgn (det3 S X R) with
  | .ok a, .ok b => .ok (if strict then decide (a > 0 ∧ b > 0) else decide (a ≥ 0 ∧ b ≥ 0))
  | .error e, _ => .error e
  | _, .error e => .error e

/-- One end of the intersection of two intervals: `X2` if it lies in `I1`,
else `X1` if it lies in `I2`. -/
def pickEnd (S : Vec) (I1 I2 : Vec × Vec) (X1 X2 : Vec) : KM (Option Vec) :=
  match inside S I1.1 I1.2 X2 false with
  | .error e => .error e
  | .ok true => .ok (some X2)
  | .ok false =>
    match inside S I2.1 I2.2 X1 false with
    | .error e => .error e
    | .ok true => .ok (some X1)
    | .ok false => .ok none

def intersect (S : Vec) (I1 I2 : Vec × Vec) : KM (Option (Vec × Vec)) :=
  match pickEnd S I1 I2 I1.1 I2.1, pickEnd S I1 I2 I1.2 I2.2 with
  | .ok (some L), .ok (some R) =>
    match sgn (det3 S L R) with
    | .ok s => .ok (if s > 0 then some (L, R) else none)
    | .error e => .error e
  | .error e, _ => .error e
  | _, .error e => .error e
  | _, _ => .ok none

/-- Record the strip ending at its new vertex `r` as a bead, if it is one. -/
def recordBead (G : Gram) (piBeads : Bool) (ts2 : Array Nat) (pos2 : Array (Array Vec))
    (steps2 : Array Nat) (start r : Nat) (S E : Vec) (iv : Vec × Vec) (par : Int)
    (out : Array Bead) : KM (Array Bead) :=
  if par < 0 then
    if piBeads then
      match Bead.make ts2 start r G pos2 steps2 none with
      | .ok b => .ok (out.push b)
      | .error e => .error e
    else .ok out
  else if par = 0 then
    match inside S iv.1 iv.2 E true with
    | .error e => .error e
    | .ok false => .ok out
    | .ok true =>
      match Bead.make ts2 start r G pos2 steps2 none with
      | .ok b => .ok (out.push b)
      | .error e => .error e
  else .ok out

/-- Continue across the edge opposite vertex `dn` of the triangle just
entered (new vertex `r` at `E`), if some direction still crosses every edge. -/
def growBranch (G : Gram) (rec : Nat → Vec × Vec → Array Bead → KM (Array Bead))
    (S E : Vec) (P2 : Array Vec) (r dn : Nat) (iv : Vec × Vec) (out : Array Bead) :
    KM (Array Bead) :=
  match parallel G S P2[3 - r - dn]! with
  | .error e => .error e
  | .ok p =>
    if p ≠ 0 then .ok out
    else
      match visible S E P2[3 - r - dn]! with
      | .error e => .error e
      | .ok none => .ok out
      | .ok (some iv2) =>
        match intersect S iv iv2 with
        | .error e => .error e
        | .ok none => .ok out
        | .ok (some new) => rec dn new out

/-- Cross the edge opposite vertex `d` of the strip's last triangle (`_grow`). -/
def grow (L : Link) (G : Gram) (piBeads : Bool) :
    Nat → Array Nat → Array (Array Vec) → Array Nat → Nat → Vec → Nat → Vec × Vec →
      Array Bead → KM (Array Bead)
  | 0, _, _, _, _, _, _, _, _ => .error "bead strip longer than 4x the link: growth not terminating?"
  | fuel + 1, ts, pos, steps, start, S, d, iv, out =>
    let t := ts.back!
    let P2 := (step L pos.back! t d).2
    let ts2 := ts.push ((L.adj[t]!)[d]!)
    let pos2 := pos.push P2
    let steps2 := steps.push d
    let r := ((L.glue[t]!)[d]!)[d]!
    let E := P2[r]!
    match parallel G S E with
    | .error e => .error e
    | .ok par =>
      match recordBead G piBeads ts2 pos2 steps2 start r S E iv par out with
      | .error e => .error e
      | .ok out1 =>
        if par ≠ 0 then .ok out1
        else
          let rec' := fun dn new o => grow L G piBeads fuel ts2 pos2 steps2 start S dn new o
          match growBranch G rec' S E P2 r ((r + 1) % 3) iv out1 with
          | .error e => .error e
          | .ok out2 => growBranch G rec' S E P2 r ((r + 2) % 3) iv out2

/-- The beads starting at corner `d0` of link triangle `t`: its two edges, and
the strips grown across the opposite edge. -/
def enumStart (L : Link) (piBeads : Bool) (t d0 : Nat) (out : Array Bead) : KM (Array Bead) :=
  let G := L.G[t]!
  let P := E3
  let a := (d0 + 1) % 3
  let b := (d0 + 2) % 3
  let S := P[d0]!
  match Bead.make #[t] d0 a G #[P] #[] none with
  | .error e => .error e
  | .ok ba =>
    match Bead.make #[t] d0 b G #[P] #[] none with
    | .error e => .error e
    | .ok bb =>
      let out := (out.push ba).push bb
      match parallel G S P[a]!, parallel G S P[b]! with
      | .ok pa, .ok pb =>
        if pa ≠ 0 ∨ pb ≠ 0 then .ok out
        else
          match visible S P[a]! P[b]! with
          | .error e => .error e
          | .ok none => .ok out
          | .ok (some iv) => grow L G piBeads (4 * L.n + 1) #[t] #[P] #[] d0 S d0 iv out
      | .error e, _ => .error e
      | _, .error e => .error e

def enumList (L : Link) (piBeads : Bool) : List (Nat × Nat) → Array Bead → KM (Array Bead)
  | [], out => .ok out
  | p :: rest, out =>
    match enumStart L piBeads p.1 p.2 out with
    | .error e => .error e
    | .ok out' => enumList L piBeads rest out'

/-- The corners `(t, d0)` of the link triangles, in enumeration order. -/
def corners (L : Link) : List (Nat × Nat) :=
  (List.range L.n).flatMap fun t => (List.range 3).map fun d0 => (t, d0)

def enumerateBeads (L : Link) (piBeads : Bool) : KM (Array Bead) :=
  enumList L piBeads (corners L) #[]

/-! ### Junctions -/

lemma angleAdd_ok {st st' : AngleSum} {ab : MQ × MQ} (h : angleAdd st ab = .ok st') :
    st.add ab.1 ab.2 = some st' := by
  unfold angleAdd at h
  split at h
  · cases h; assumption
  · cases h

lemma angleLeq_ok {x y : AngleSum} {b : Bool} (h : angleLeq x y = .ok b) : x.leq y = some b := by
  unfold angleLeq at h
  split at h
  · cases h; assumption
  · cases h

/-- The fan walk around the junction's link vertex: from corner `c` of
triangle `t`, crossing the edge opposite `ent`, through the corners of the
triangles at that link vertex until the outgoing bead's corner `(t0, d0)`.
Returns the exact sum and the list of angles added (stopping as soon as the
sum is `≥ π`, as `kernel.junction_side` does). -/
def fanWalk (L : Link) (t0 d0 : Nat) (outG : Gram) (outDir : V3) :
    Nat → Nat → Nat → Nat → Nat → AngleSum → List (MQ × MQ) → KM (AngleSum × List (MQ × MQ))
  | 0, _, _, _, _, _, _ => .error "junction fan walk did not reach the outgoing triangle"
  | fuel + 1, t, c, ent, ex, acc, l =>
    if atLeast acc 1 then .ok (acc, l)
    else
      let g := (L.glue[t]!)[ent]!
      let t' := (L.adj[t]!)[ent]!
      let c' := g[c]!
      let ent' := g[ex]!
      let ex' := 3 - c' - ent'
      if (t', c') = (t0, d0) then
        match angleAdd acc (tangentAngle outG E3[d0]! E3[ent']! outDir) with
        | .ok acc' => .ok (acc', l ++ [tangentAngle outG E3[d0]! E3[ent']! outDir])
        | .error e => .error e
      else
        match angleAdd acc (tangentAngle L.G[t']! E3[c']! E3[ent']! E3[ex']!) with
        | .ok acc' => fanWalk L t0 d0 outG outDir fuel t' c' ent' ex' acc'
            (l ++ [tangentAngle L.G[t']! E3[c']! E3[ent']! E3[ex']!])
        | .error e => .error e

lemma sum_append_angle (l : List (MQ × MQ)) (ab : MQ × MQ) :
    ((l ++ [ab]).map fun p => angleOf p.1 p.2).sum =
      (l.map fun p => angleOf p.1 p.2).sum + angleOf ab.1 ab.2 := by
  simp

theorem fanWalk_ok {L : Link} {t0 d0 : Nat} {outG : Gram} {outDir : V3} :
    ∀ (fuel t c ent ex : Nat) (acc : AngleSum) (l : List (MQ × MQ)) (r : AngleSum × List (MQ × MQ)),
      acc.Rep ((l.map fun p => angleOf p.1 p.2).sum) →
      fanWalk L t0 d0 outG outDir fuel t c ent ex acc l = .ok r →
      r.1.Rep ((r.2.map fun p => angleOf p.1 p.2).sum) := by
  intro fuel
  induction fuel with
  | zero => intro t c ent ex acc l r _ h; simp [fanWalk] at h
  | succ fuel ih =>
    intro t c ent ex acc l r hr h
    unfold fanWalk at h
    dsimp only at h
    split_ifs at h with hat hcorner
    · cases h; exact hr
    · split at h
      · rename_i acc' hadd
        cases h
        rw [sum_append_angle]
        exact AngleSum.add_sound hr (angleAdd_ok hadd)
      · cases h
    · split at h
      · rename_i acc' hadd
        refine ih _ _ _ _ acc' _ r ?_ h
        rw [sum_append_angle]
        exact AngleSum.add_sound hr (angleAdd_ok hadd)
      · cases h

/-- Is the outgoing direction `X` strictly inside the corner at `W` between
`u` and the ray toward `Py`? (Their angles add up exactly.) -/
def insideCorner (G : Gram) (W u X Py : V3) : KM Bool :=
  match angleSum #[tangentAngle G W u X, tangentAngle G W X Py],
      angleSum #[tangentAngle G W u Py] with
  | .ok s12, .ok s3 =>
    if s12.h = s3.h then
      match twSgn (Tw.sub (Tw.mul s12.im s3.re) (Tw.mul s12.re s3.im)) with
      | .ok 0 =>
        match twSgn (Tw.add (Tw.mul s12.re s3.re) (Tw.mul s12.im s3.im)) with
        | .ok re => .ok (decide (re > 0))
        | .error e => .error e
      | .ok _ => .ok false
      | .error e => .error e
    else .ok false
  | .error e, _ => .error e
  | _, .error e => .error e

/-- The angle at the junction `bin → bout` on side `k ∈ {1, 2}`, with the
corner angles it sums. -/
def junctionFan (L : Link) (bin bout : Bead) (k : Nat) : KM (AngleSum × List (MQ × MQ)) :=
  let tm := bin.ts.back!
  let rm := bin.stop
  let t0 := bout.ts[0]!
  let d0 := bout.start
  let P := bin.pos.back!
  let W := P[rm]!
  let u := bin.inDir
  let d := (rm + k) % 3
  let y := 3 - rm - d
  let first := tangentAngle bin.G W u P[y]!
  let general : KM (AngleSum × List (MQ × MQ)) :=
    match angleSum #[first] with
    | .ok acc => fanWalk L t0 d0 bout.G bout.outDir (4 * L.n + 2) tm rm d y acc [first]
    | .error e => .error e
  if (tm, rm) = (t0, d0) then
    -- the outgoing direction, in this triangle's developed frame
    let X : V3 := V3.add (V3.add (V3.scale bout.outDir.x P[0]!) (V3.scale bout.outDir.y P[1]!))
      (V3.scale bout.outDir.z P[2]!)
    match insideCorner bin.G W u X P[y]! with
    | .ok true =>
      match angleSum #[tangentAngle bin.G W u X] with
      | .ok s => .ok (s, [tangentAngle bin.G W u X])
      | .error e => .error e
    | .ok false => general
    | .error e => .error e
  else general

/-- The exact angle at a junction side is the sum of the corner angles listed. -/
theorem junctionFan_ok {L : Link} {bin bout : Bead} {k : Nat} {r : AngleSum × List (MQ × MQ)}
    (h : junctionFan L bin bout k = .ok r) :
    r.1.Rep ((r.2.map fun p => angleOf p.1 p.2).sum) := by
  have hgen : ∀ (first : MQ × MQ) (t0 d0 tm rm d y : Nat),
      (match angleSum #[first] with
        | .ok acc => fanWalk L t0 d0 bout.G bout.outDir (4 * L.n + 2) tm rm d y acc [first]
        | .error e => .error e) = .ok r →
      r.1.Rep ((r.2.map fun p => angleOf p.1 p.2).sum) := by
    intro first t0 d0 tm rm d y hg
    split at hg
    · rename_i acc hacc
      have := angleSum_ok hacc
      exact fanWalk_ok _ _ _ _ _ acc [first] r (by simpa using this) hg
    · cases hg
  unfold junctionFan at h
  dsimp only at h
  split_ifs at h with hsame
  · split at h
    · split at h
      · rename_i s hs
        cases h
        simpa using angleSum_ok hs
      · cases h
    · exact hgen _ _ _ _ _ _ _ h
    · cases h
  · exact hgen _ _ _ _ _ _ _ h

/-- Every angle the fan walk adds is a corner angle at the junction's link
vertex: the whole corner of a link triangle there (between its two edges at
that vertex), or, in the outgoing triangle, the part from its entry edge to
the outgoing direction. -/
theorem fanWalk_corners {L : Link} (hg : glueOk L = true) {t0 d0 : Nat} {outG : Gram}
    {outDir : V3} :
    ∀ (fuel t c ent ex : Nat) (acc : AngleSum) (l : List (MQ × MQ)) (r : AngleSum × List (MQ × MQ)),
      t < L.n → Corner3 c ent ex →
      fanWalk L t0 d0 outG outDir fuel t c ent ex acc l = .ok r →
      ∃ extra, r.2 = l ++ extra ∧ ∀ p ∈ extra, ∃ t' c' a b, t' < L.n ∧ Corner3 c' a b ∧
        (L.lv[t']!)[c']! = (L.lv[t]!)[c]! ∧
        (p = tangentAngle L.G[t']! E3[c']! E3[a]! E3[b]! ∨
          ((t', c') = (t0, d0) ∧ p = tangentAngle outG E3[d0]! E3[a]! outDir)) := by
  intro fuel
  induction fuel with
  | zero => intro t c ent ex acc l r _ _ h; simp [fanWalk] at h
  | succ fuel ih =>
    intro t c ent ex acc l r ht hc h
    obtain ⟨hn, hc', hW, -, -⟩ := fan_step hg ht hc
    simp only [fanNext] at hn hc' hW
    unfold fanWalk at h
    dsimp only at h
    split_ifs at h with hat hcorner
    · cases h; exact ⟨[], by simp, by simp⟩
    · split at h
      · rename_i acc' _
        cases h
        refine ⟨[_], rfl, ?_⟩
        intro p hp
        simp only [List.mem_singleton] at hp
        subst hp
        refine ⟨_, _, _, 3 - ((L.glue[t]!)[ent]!)[c]! - ((L.glue[t]!)[ent]!)[ex]!, hn, hc', hW, ?_⟩
        right
        exact ⟨hcorner, rfl⟩
      · cases h
    · split at h
      · rename_i acc' _
        obtain ⟨extra, hr, hall⟩ := ih _ _ _ _ acc' _ r hn hc' h
        refine ⟨_ :: extra, hr.trans (List.append_assoc _ _ _), ?_⟩
        intro p hp
        rcases List.mem_cons.mp hp with rfl | hp
        · exact ⟨_, _, _, _, hn, hc', hW, Or.inl rfl⟩
        · obtain ⟨t', c', a, b, h1, h2, h3, h4⟩ := hall p hp
          exact ⟨t', c', a, b, h1, h2, h3.trans hW, h4⟩
      · cases h

def junctionSide (L : Link) (bin bout : Bead) (k : Nat) : KM AngleSum :=
  match junctionFan L bin bout k with
  | .ok r => .ok r.1
  | .error e => .error e

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

/-- The strip never goes straight back across the edge it just crossed. -/
def noBacktrack (L : Link) (spec : BeadSpec) : Bool :=
  (List.range spec.steps.size).all fun k =>
    k == 0 || spec.steps[k]! != ((L.glue[spec.ts[k - 1]!]!)[spec.steps[k - 1]!]!)[spec.steps[k - 1]!]!

lemma noBacktrack_spec {L : Link} {spec : BeadSpec} (h : noBacktrack L spec = true) (k : Nat)
    (hk0 : 0 < k) (hk : k < spec.steps.size) :
    spec.steps[k]! ≠ ((L.glue[spec.ts[k - 1]!]!)[spec.steps[k - 1]!]!)[spec.steps[k - 1]!]! := by
  unfold noBacktrack at h
  simp only [List.all_eq_true, List.mem_range, Bool.or_eq_true, beq_iff_eq, bne_iff_ne,
    ne_eq] at h
  exact (h k hk).resolve_left (by omega)

def checkBead (L : Link) (spec : BeadSpec) : KM Bead :=
  guardK (!(spec.ts.size = 0 || spec.ts.any (· ≥ L.n) || spec.start > 2 || spec.stop > 2))
    "malformed bead" >>= fun _ =>
  develop L spec.ts spec.steps >>= fun pos =>
  guardK (noBacktrack L spec) "bead strip turns back" >>= fun _ =>
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
    b.ts = spec.ts ∧ b.steps = spec.steps ∧ b.start = spec.start ∧ b.stop = spec.stop ∧
      (∀ t ∈ spec.ts.toList, t < L.n) ∧ spec.ts.size ≠ 0 ∧ spec.start < 3 ∧
      develop L spec.ts spec.steps = .ok b.pos ∧ noBacktrack L spec = true ∧
      b.S = (b.pos[0]!)[b.start]! ∧ b.E = b.pos.back![b.stop]! ∧
      (b.pi = false → b.outDir = b.E) ∧
      ∀ k ∈ List.range (b.ts.size - 1),
        crosses b.G b.S b.outDir ((b.pos[k]!)[(others3 b.steps[k]!).1]!)
          ((b.pos[k]!)[(others3 b.steps[k]!).2]!) b.pi = .ok true := by
  unfold checkBead at h
  obtain ⟨_, hm, h⟩ := bind_ok h
  obtain ⟨pos, hdev, h⟩ := bind_ok h
  obtain ⟨_, hnb, h⟩ := bind_ok h
  obtain ⟨par, -, h⟩ := bind_ok h
  obtain ⟨_, -, h⟩ := bind_ok h
  obtain ⟨Q, -, h⟩ := bind_ok h
  obtain ⟨b', hb, h⟩ := bind_ok h
  obtain ⟨_, hc, h⟩ := bind_ok h
  cases h
  obtain ⟨hts, hG, hposb, hst, hs, he, hS, hE, hout⟩ := Bead.make_ok hb
  have hm' := guardK_ok hm
  simp only [Bool.not_eq_true', Bool.or_eq_false_iff, decide_eq_false_iff_not, not_lt,
    Array.any_eq_false, not_le] at hm'
  refine ⟨hts, hst, hs, he, ?_, ?_, ?_, by rw [hposb]; exact hdev, guardK_ok hnb,
    by rw [hS, hposb, hs], by rw [hE, hposb, he], ?_, ?_⟩
  · intro t ht
    obtain ⟨i, hi, rfl⟩ := List.mem_iff_getElem.mp ht
    simpa using hm'.1.1.2 i (by simpa using hi)
  · exact hm'.1.1.1
  · omega
  · intro hpi; rw [hout, hpi]; rfl
  · have := checkCrossings_ok hc
    rw [hts, hG, hposb, hst, hS]
    exact this


/-! ### Beads cross their strip in order -/

lemma others3_pair {d p q : Nat} (hd : d < 3) (hp : p < 3) (hq : q < 3) (hpd : p ≠ d)
    (hqd : q ≠ d) (hpq : p ≠ q) : others3 d = (p, q) ∨ others3 d = (q, p) := by
  unfold others3; split_ifs <;> simp only [Prod.mk.injEq] <;> omega

lemma others3_mem {d i : Nat} (hd : d < 3) (hi : i < 3) (hid : i ≠ d) :
    i = (others3 d).1 ∨ i = (others3 d).2 := by
  unfold others3; split_ifs <;> simp <;> omega

/-- The vertex `stepWith` adds. -/
def stepNew (tr : Trans) (P : Array V3) (d : Nat) : V3 :=
  V3.sub (V3.add (V3.scaleQ tr.q.1 P[(others3 d).1]!) (V3.scaleQ tr.q.2 P[(others3 d).2]!))
    (V3.scale tr.mu (V3.sub (V3.sub P[d]! (V3.scaleQ tr.p.1 P[(others3 d).1]!))
      (V3.scaleQ tr.p.2 P[(others3 d).2]!)))

/-- Where `stepWith` puts things: the shared edge's ends stay, the new vertex
is `stepNew`. -/
lemma stepWith_idx (tr : Trans) (P : Array V3) {d : Nat} (hd : d < 3)
    (hperm : isPerm3 tr.g = true) :
    (stepWith tr P d)[tr.g[(others3 d).1]!]! = P[(others3 d).1]! ∧
      (stepWith tr P d)[tr.g[(others3 d).2]!]! = P[(others3 d).2]! ∧
      (stepWith tr P d)[tr.g[d]!]! = stepNew tr P d := by
  obtain ⟨ha3, hb3, hab, had, hbd⟩ := others3_spec hd
  obtain ⟨gab3, gab⟩ := isPerm3_spec hperm ha3 hb3 hab
  obtain ⟨-, gad⟩ := isPerm3_spec hperm ha3 hd had
  obtain ⟨gb3, gbd⟩ := isPerm3_spec hperm hb3 hd hbd
  obtain ⟨gd3, -⟩ := isPerm3_spec hperm hd ha3 (Ne.symm had)
  have hE : ∀ i < 3, (stepWith tr P d)[i]! =
      if i = tr.g[(others3 d).1]! then P[(others3 d).1]! else
        if i = tr.g[(others3 d).2]! then P[(others3 d).2]! else stepNew tr P d :=
    fun i hi => by interval_cases i <;> rfl
  exact ⟨by rw [hE _ gab3, if_pos rfl], by rw [hE _ gb3, if_neg (Ne.symm gab), if_pos rfl],
    by rw [hE _ gd3, if_neg (Ne.symm gad), if_neg (Ne.symm gbd)]⟩

/-- The new vertex lies across the shared edge's plane from the old one. -/
lemma stepNew_det (tr : Trans) (P : Array V3) (d : Nat) :
    rdet (P[(others3 d).1]!).val (P[(others3 d).2]!).val (stepNew tr P d).val =
      -tr.mu.val * rdet (P[(others3 d).1]!).val (P[(others3 d).2]!).val (P[d]!).val := by
  simp only [stepNew, V3.val_sub, V3.val_add, V3.val_scaleQ, V3.val_scale]
  simp only [rdet_eq, Pi.sub_apply, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  ring

lemma mkTrans_tr {G : Array Gram} {adj : Array (Array Nat)} {glue : Array (Array (Array Nat))}
    {t d : Nat} {tr : Trans} (h : mkTrans G adj glue t d = .ok tr) :
    tr.g = (glue[t]!)[d]! ∧ isPerm3 tr.g = true ∧ 0 < tr.mu.val := by
  unfold mkTrans at h
  simp only at h
  split_ifs at h with h1 h2 h3 <;> try (cases h; done)
  cases h
  simp only [Bool.not_eq_true', Bool.and_eq_false_iff, Bool.not_eq_false, not_or] at h1
  push Not at h3
  obtain ⟨-, hh, hh2⟩ := h3
  refine ⟨rfl, h1.1.1.1.1, ?_⟩
  simp only
  rw [MQ.val_sqrtRat _ (le_of_lt (div_pos hh2 hh))]
  exact Real.sqrt_pos.mpr (by exact_mod_cast div_pos hh2 hh)

/-- `develop` is a sequence of `step`s along the strip. -/
lemma developFrom_step {L : Link} :
    ∀ (ds : List Nat) (t : Nat) (ts : List Nat) (P : Array V3) (out : List (Array V3)),
      developFrom L P (t :: ts) ds = .ok out →
      out.length = ds.length ∧ ts.length = ds.length ∧
      ∀ k < ds.length, (P :: out)[k + 1]! =
          stepWith (L.tr (t :: ts)[k]! ds[k]!) (P :: out)[k]! ds[k]! ∧
        ds[k]! < 3 ∧ (L.adj[(t :: ts)[k]!]!)[ds[k]!]! = (t :: ts)[k + 1]! := by
  intro ds
  induction ds with
  | nil =>
    intro t ts P out h
    cases ts with
    | nil => simp [developFrom] at h; subst h; simp
    | cons t' ts => simp [developFrom] at h
  | cons d ds ih =>
    intro t ts P out h
    cases ts with
    | nil => simp [developFrom] at h
    | cons t' ts =>
      unfold developFrom at h
      split_ifs at h with hbad
      push Not at hbad
      split at h
      · rename_i rest hrest
        cases h
        obtain ⟨hl1, hl2, hk⟩ := ih t' ts _ rest hrest
        refine ⟨by simp [hl1], by simp [hl2], ?_⟩
        intro k hk'
        cases k with
        | zero => simp [step, hbad.2]; omega
        | succ k =>
          have := hk k (by simpa using hk')
          simpa using this
      · cases h

lemma develop_step {L : Link} {ts steps : Array Nat} {pos : Array (Array V3)}
    (h : develop L ts steps = .ok pos) :
    pos.size = ts.size ∧ steps.size + 1 = ts.size ∧ ∀ k < steps.size,
      pos[k + 1]! = stepWith (L.tr ts[k]! steps[k]!) pos[k]! steps[k]! ∧ steps[k]! < 3 ∧
        (L.adj[ts[k]!]!)[steps[k]!]! = ts[k + 1]! := by
  unfold develop at h
  split_ifs at h with hsz
  split at h
  · rename_i rest hrest
    cases h
    obtain ⟨t0, tl, htl⟩ : ∃ t0 tl, ts.toList = t0 :: tl := by
      cases hc : ts.toList with
      | nil => have : ts.size = 0 := by simpa using congrArg List.length hc
               omega
      | cons t0 tl => exact ⟨t0, tl, rfl⟩
    rw [htl] at hrest
    obtain ⟨hl1, hl2, hk⟩ := developFrom_step _ _ _ _ _ hrest
    have hts : ∀ i : Nat, ts[i]! = (t0 :: tl)[i]! := by
      intro i; rw [← htl, Array.getElem!_toList]
    have hst : ∀ i : Nat, steps[i]! = steps.toList[i]! := by intro i; rw [Array.getElem!_toList]
    refine ⟨by simp [hl1]; omega, by omega, ?_⟩
    intro k hk'
    have := hk k (by simpa using hk')
    rw [hts, hts, hst]
    simpa using this
  · cases h

/-- An accepted bead (not a π-bead) is a geodesic through its strip in order:
the segment from `S` toward `E` meets the open edge `k` of the strip at a
parameter `τ k > 0`, and these parameters increase along the strip. -/
theorem checkBead_order {geo : Geometry} {iv : Nat} {L : Link} (hL : Link.make geo iv = .ok L)
    {spec : BeadSpec} {b : Bead} (h : checkBead L spec = .ok b) (hpi : b.pi = false) :
    ∃ τ : ℕ → ℝ,
      (∀ k, k + 1 < b.ts.size → 0 < τ k ∧ ∃ a1 b1 : ℝ, 0 < a1 ∧ 0 < b1 ∧
        b.S.val + τ k • b.E.val =
          a1 • ((b.pos[k]!)[(others3 b.steps[k]!).1]!).val +
            b1 • ((b.pos[k]!)[(others3 b.steps[k]!).2]!).val) ∧
      (∀ k, k + 2 < b.ts.size → τ k < τ (k + 1)) := by
  obtain ⟨hts, hst, hs, he, hn, hsz, hstart, hdev, hnb, hS, hE, hout, hcr⟩ := checkBead_ok h
  rw [← hts, ← hst] at hdev
  obtain ⟨hpsz, hsteps, hdstep⟩ := develop_step hdev
  have hout' : b.outDir = b.E := hout hpi
  let P : ℕ → Array V3 := fun k => b.pos[k]!
  let d : ℕ → ℕ := fun k => b.steps[k]!
  let A : ℕ → (Fin 3 → ℝ) := fun k => ((P k)[(others3 (d k)).1]!).val
  let B : ℕ → (Fin 3 → ℝ) := fun k => ((P k)[(others3 (d k)).2]!).val
  let D : ℕ → (Fin 3 → ℝ) := fun k => ((P k)[d k]!).val
  -- crossing data for every interior edge
  have data : ∀ k, ∃ τ a1 b1 : ℝ, k + 1 < b.ts.size → 0 < τ ∧ 0 < a1 ∧ 0 < b1 ∧
      b.S.val + τ • b.E.val = a1 • A k + b1 • B k ∧ rdet (A k) (B k) b.S.val ≠ 0 := by
    intro k
    by_cases hk : k + 1 < b.ts.size
    · have hc := hcr k (List.mem_range.mpr (by omega))
      rw [hout', hpi] at hc
      obtain ⟨α, β, γ, δ, hα, hβ, hγ, hδ, hcomb, hdS⟩ := crosses_ok_minor hc
      refine ⟨δ / γ, α / γ, β / γ, fun _ => ⟨div_pos hδ hγ, div_pos hα hγ, div_pos hβ hγ, ?_, hdS⟩⟩
      ext i
      have hi := congrFun hcomb i
      simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul] at hi ⊢
      field_simp
      linarith
    · exact ⟨0, 0, 0, fun h => absurd h hk⟩
  choose τ a1 b1 hdat using data
  -- one step along the strip
  have step_k : ∀ k, k + 2 < b.ts.size →
      0 < rdet (A k) (B k) b.S.val * rdet (A k) (B k) (D k) →
      τ k < τ (k + 1) ∧
        0 < rdet (A (k + 1)) (B (k + 1)) b.S.val * rdet (A (k + 1)) (B (k + 1)) (D (k + 1)) := by
    intro k hk hside
    have hks : k < b.steps.size := by omega
    have hks1 : k + 1 < b.steps.size := by omega
    obtain ⟨hPk, hdk, -⟩ := hdstep k hks
    obtain ⟨-, hdk1, -⟩ := hdstep (k + 1) hks1
    have htk : b.ts[k]! < L.n := by
      apply hn
      have hk' : k < b.ts.size := by omega
      rw [← hts, getElem!_pos b.ts k hk']
      simp
    have hm := Link.make_trans hL _ htk _ hdk
    obtain ⟨hg, hperm, hμ⟩ := mkTrans_tr hm
    obtain ⟨iA, iB, iN⟩ := stepWith_idx (L.tr b.ts[k]! b.steps[k]!) (P k) hdk hperm
    have hdet := stepNew_det (L.tr b.ts[k]! b.steps[k]!) (P k) (d k)
    obtain ⟨ha3, hb3, hab, had, hbd⟩ := others3_spec hdk
    obtain ⟨gab3, gab⟩ := isPerm3_spec hperm ha3 hb3 hab
    obtain ⟨-, gad⟩ := isPerm3_spec hperm ha3 hdk had
    obtain ⟨gb3, gbd⟩ := isPerm3_spec hperm hb3 hdk hbd
    obtain ⟨gd3, -⟩ := isPerm3_spec hperm hdk ha3 (Ne.symm had)
    set g := (L.tr b.ts[k]! b.steps[k]!).g with hgdef
    -- no backtracking: the next step does not cross back
    have hnb' : b.steps[k + 1]! ≠ g[b.steps[k]!]! := by
      have := noBacktrack_spec hnb (k + 1) (by omega) (by rw [← hst]; omega)
      simp only [Nat.add_sub_cancel] at this
      rw [hg]
      simpa [hts, hst] using this
    have hP1 : P (k + 1) = stepWith (L.tr b.ts[k]! b.steps[k]!) (P k) (d k) := hPk
    obtain ⟨hτk, ha1k, hb1k, hck, -⟩ := hdat k (by omega)
    obtain ⟨hτk1, ha2, hb2, hck1, -⟩ := hdat (k + 1) hk
    set N := (stepNew (L.tr b.ts[k]! b.steps[k]!) (P k) (d k)).val
    have vA : ((P (k + 1))[g[(others3 b.steps[k]!).1]!]!).val = A k := by rw [hP1, iA]
    have vB : ((P (k + 1))[g[(others3 b.steps[k]!).2]!]!).val = B k := by rw [hP1, iB]
    have vN : ((P (k + 1))[g[b.steps[k]!]!]!).val = N := by rw [hP1, iN]
    have hcase : (A (k + 1) = A k ∧ B (k + 1) = N ∧ D (k + 1) = B k) ∨
        (A (k + 1) = N ∧ B (k + 1) = A k ∧ D (k + 1) = B k) ∨
        (A (k + 1) = B k ∧ B (k + 1) = N ∧ D (k + 1) = A k) ∨
        (A (k + 1) = N ∧ B (k + 1) = B k ∧ D (k + 1) = A k) := by
      have hcov : b.steps[k + 1]! = g[(others3 b.steps[k]!).2]! ∨
          b.steps[k + 1]! = g[(others3 b.steps[k]!).1]! := by
        omega
      rcases hcov with hdB | hdA
      · rcases others3_pair hdk1 gab3 gd3 (by omega) (Ne.symm hnb') gad with hp | hp
        · left
          refine ⟨?_, ?_, ?_⟩
          · show ((P (k + 1))[(others3 b.steps[k + 1]!).1]!).val = A k; rw [hp]; exact vA
          · show ((P (k + 1))[(others3 b.steps[k + 1]!).2]!).val = N; rw [hp]; exact vN
          · show ((P (k + 1))[b.steps[k + 1]!]!).val = B k; rw [hdB]; exact vB
        · right; left
          refine ⟨?_, ?_, ?_⟩
          · show ((P (k + 1))[(others3 b.steps[k + 1]!).1]!).val = N; rw [hp]; exact vN
          · show ((P (k + 1))[(others3 b.steps[k + 1]!).2]!).val = A k; rw [hp]; exact vA
          · show ((P (k + 1))[b.steps[k + 1]!]!).val = B k; rw [hdB]; exact vB
      · rcases others3_pair hdk1 gb3 gd3 (by omega) (Ne.symm hnb') gbd with hp | hp
        · right; right; left
          refine ⟨?_, ?_, ?_⟩
          · show ((P (k + 1))[(others3 b.steps[k + 1]!).1]!).val = B k; rw [hp]; exact vB
          · show ((P (k + 1))[(others3 b.steps[k + 1]!).2]!).val = N; rw [hp]; exact vN
          · show ((P (k + 1))[b.steps[k + 1]!]!).val = A k; rw [hdA]; exact vA
        · right; right; right
          refine ⟨?_, ?_, ?_⟩
          · show ((P (k + 1))[(others3 b.steps[k + 1]!).1]!).val = N; rw [hp]; exact vN
          · show ((P (k + 1))[(others3 b.steps[k + 1]!).2]!).val = B k; rw [hp]; exact vB
          · show ((P (k + 1))[b.steps[k + 1]!]!).val = A k; rw [hdA]; exact vA
    exact order_next hμ hdet hside hτk ha1k hb1k hck hcase ha2 hb2 hck1
  -- the start point is on the right side of every edge
  have side : ∀ k, k + 1 < b.ts.size → 0 < rdet (A k) (B k) b.S.val * rdet (A k) (B k) (D k) := by
    intro k
    induction k with
    | zero =>
      intro h1
      obtain ⟨-, -, -, -, hdS⟩ := hdat 0 h1
      have hd0 : d 0 < 3 := (hdstep 0 (by omega)).2.1
      have hstart' : b.start < 3 := by rw [hs]; exact hstart
      have hS0 : b.S.val = ((P 0)[b.start]!).val := by rw [hS]
      have hdst : d 0 = b.start := by
        by_contra hne
        rcases others3_mem hd0 hstart' (Ne.symm hne) with hm | hm
        · apply hdS
          show rdet ((P 0)[(others3 (d 0)).1]!).val ((P 0)[(others3 (d 0)).2]!).val b.S.val = 0
          rw [hS0, hm]; exact rdet_self_left _ _
        · apply hdS
          show rdet ((P 0)[(others3 (d 0)).1]!).val ((P 0)[(others3 (d 0)).2]!).val b.S.val = 0
          rw [hS0, hm]; exact rdet_self_right _ _
      have hD0 : D 0 = b.S.val := by show ((P 0)[d 0]!).val = b.S.val; rw [hdst, hS0]
      rw [hD0]
      exact mul_self_pos.mpr hdS
    | succ k ih =>
      intro h1
      exact (step_k k (by omega) (ih (by omega))).2
  refine ⟨τ, fun k hk => ?_, fun k hk => (step_k k hk (side k (by omega))).1⟩
  obtain ⟨hτk, ha, hb, hc, -⟩ := hdat k hk
  exact ⟨hτk, a1 k, b1 k, ha, hb, hc⟩

inductive Reason where
  | invalid (tet : Nat)
  | edge (e : Nat)
  | vertex (v : Nat) (necklace : Array BeadSpec)

/-- A certificate that the metric is not locally CAT(0) (`vlink-certificate/2`). -/
structure FalseCert where
  isosig : String
  sqlens : Array ℚ
  reason : Reason

def startV (L : Link) (b : Bead) : Nat := (L.lv[b.ts[0]!]!)[b.start]!
def endV (L : Link) (b : Bead) : Nat := (L.lv[b.ts.back!]!)[b.stop]!

/-- Every tetrahedron is valid. -/
def tetsValidK (T : Tri) (geo : Geometry) (msg : Nat → String) : KM Unit :=
  forallK (List.range T.n) fun t => if geo.tetValid t then .ok () else .error (msg t)

lemma tetsValidK_ok {T : Tri} {geo : Geometry} {msg : Nat → String}
    (h : tetsValidK T geo msg = .ok ()) : ∀ t < T.n, geo.tetValid t = true := by
  intro t ht
  have := forallK_ok h t (List.mem_range.mpr ht)
  by_contra hv
  rw [if_neg hv] at this
  cases this

/-- Both sides of the junction `b1 → b2` have angle `≥ π`. -/
def junctionPasses (L : Link) (b1 b2 : Bead) (k : Nat) : KM Unit :=
  forallK [1, 2] fun side =>
    match junctionSide L b1 b2 side with
    | .ok acc => if atLeast acc 1 then .ok () else .error s!"junction {k}: angle < pi on side {side}"
    | .error e => .error e

/-- The beads meet in a cycle, with every junction passing. -/
def necklaceOk (L : Link) (beads : List Bead) : KM Unit :=
  forallK (List.range beads.length) fun k =>
    if endV L beads[k]! ≠ startV L beads[(k + 1) % beads.length]! then
      .error s!"beads {k} and {(k + 1) % beads.length} do not meet"
    else junctionPasses L beads[k]! beads[(k + 1) % beads.length]! k

/-- Check a certificate that the metric is not locally CAT(0). -/
def checkFalse (c : FalseCert) : KM String :=
  match Tri.ofIsoSig c.isosig with
  | .error e => .error e
  | .ok T =>
  match Geometry.make T c.sqlens with
  | .error e => .error e
  | .ok geo =>
  match c.reason with
  | .invalid t =>
    if t ≥ T.n || geo.tetValid t then .error s!"tetrahedron {t} is valid"
    else .ok s!"tetrahedron {t} is degenerate"
  | .edge e =>
    match tetsValidK T geo (fun t => s!"tetrahedron {t} is degenerate (use reason \"invalid\")") with
    | .error err => .error err
    | .ok () =>
      if e ≥ T.edges.size then .error "bad edge"
      else match geo.edgeAngle e with
        | .error err => .error err
        | .ok acc =>
          if atLeast acc 2 then .error s!"edge {e} has dihedral sum >= 2*pi"
          else .ok s!"edge {e} has dihedral sum < 2*pi"
  | .vertex iv necklace =>
    match tetsValidK T geo (fun t => s!"tetrahedron {t} is degenerate (use reason \"invalid\")") with
    | .error err => .error err
    | .ok () =>
      if iv ≥ T.corners.size then .error "bad vertex"
      else match Link.make geo iv with
        | .error err => .error err
        | .ok L =>
          match mapK (checkBead L) necklace.toList with
          | .error err => .error err
          | .ok beads =>
            if beads.isEmpty then .error "empty necklace"
            else match totalLength beads.toArray with
              | .error err => .error err
              | .ok tot =>
                if atLeast tot 2 then .error "necklace length >= 2*pi"
                else match necklaceOk L beads with
                  | .error err => .error err
                  | .ok () =>
                    .ok s!"vertex {iv}: necklace of {beads.length} beads, length < 2*pi, all junctions >= pi"

/-! ### True certificates (docs/true_certificates.md) -/

def LEMMA_4_9 : String := "Elder-McCammond-Meier 2003, Lemma 4.9"

lemma angleOf_nonneg (a b2 : MQ) : 0 ≤ angleOf a b2 :=
  Complex.arg_nonneg_iff.mpr (Real.sqrt_nonneg _)


def junctionEntry (L : Link) (beads : Array Bead) (x y : Nat) : KM (Nat × Bool) :=
  if beads[x]!.pi || beads[y]!.pi then .ok (y, true)
  else match junctionOk L beads[x]! beads[y]! with
    | .ok ok => .ok (y, ok)
    | .error e => .error e

/-- Junctions out of bead `x`: every bead `y` starting where `x` ends, with
whether the junction passes (angle `≥ π` on both sides). Pairs involving a
π-bead count as passing (conservative). -/
def junctionRow (L : Link) (beads : Array Bead) (x : Nat) : KM (List (Nat × Bool)) :=
  mapK (junctionEntry L beads x)
    ((List.range beads.size).filter (fun y => startV L beads[y]! == endV L beads[x]!))

/-- `(y, true)` in the row of `x`: `y` starts where `x` ends and the junction passes. -/
theorem junctionRow_mem {L : Link} {beads : Array Bead} {x : Nat} {row : List (Nat × Bool)}
    (h : junctionRow L beads x = .ok row) {y : Nat} (hy : (y, true) ∈ row) :
    y < beads.size ∧ startV L beads[y]! = endV L beads[x]! ∧
      (beads[x]!.pi = true ∨ beads[y]!.pi = true ∨
        junctionOk L beads[x]! beads[y]! = .ok true) := by
  unfold junctionRow at h
  obtain ⟨a, ha, hga⟩ := forall2_mem_right (mapK_ok h) hy
  simp only [List.mem_filter, List.mem_range, beq_iff_eq] at ha
  unfold junctionEntry at hga
  split_ifs at hga with hpi
  · simp only [Except.ok.injEq, Prod.mk.injEq, and_true] at hga
    subst hga
    simp only [Bool.or_eq_true] at hpi
    exact ⟨ha.1, ha.2, by tauto⟩
  · split at hga
    · rename_i okv hj
      simp only [Except.ok.injEq, Prod.mk.injEq] at hga
      obtain ⟨rfl, rfl⟩ := hga
      exact ⟨ha.1, ha.2, Or.inr (Or.inr hj)⟩
    · cases hga

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

def beadLen (beads : Array Bead) (y : Nat) : MQ × MQ := beadLength beads[y]!

/-- One tree node `x` of the tree of `s`: its chain, checked to run along
passing junctions, and the exact sum of the chain's lengths. -/
def nodeEntry (beads : Array Bead) (Jt : List (List (Nat × Bool)))
    (tree : Std.HashMap Nat Nat) (s : Nat) (xp : Nat × Nat) : KM (Nat × AngleSum) :=
  if xp.1 ≤ s || xp.1 ≥ beads.size then .error s!"tree of {s}: bad bead index"
  else match chainOf tree s xp.1 with
    | .error e => .error e
    | .ok c =>
      match forallK (c.toList.zip c.toList.tail) (fun (uw : Nat × Nat) =>
          if (Jt.getD uw.1 []).contains (uw.2, true) then .ok ()
          else .error s!"tree of {s}: {uw.1} -> {uw.2} is not a passing junction") with
      | .error e => .error e
      | .ok () =>
        match angleSum (c.map (beadLen beads)) with
        | .ok sx => .ok (xp.1, sx)
        | .error e => .error e

lemma nodeEntry_rep {beads : Array Bead} {Jt : List (List (Nat × Bool))}
    {tree : Std.HashMap Nat Nat} {s : Nat} {xp : Nat × Nat} {r : Nat × AngleSum}
    (h : nodeEntry beads Jt tree s xp = .ok r) : ∃ Θ, r.2.Rep Θ := by
  unfold nodeEntry at h
  split_ifs at h
  split at h
  · cases h
  · split at h
    · cases h
    · split at h
      · rename_i sx hs
        cases h
        exact ⟨_, angleSum_ok hs⟩
      · cases h

/-- Tree nodes of the tree rooted at `s`, with the exact length sums of their
chains (`s` itself first, with the length of bead `s`). -/
def nodeSums (beads : Array Bead) (Jt : List (List (Nat × Bool)))
    (tree : Std.HashMap Nat Nat) (s : Nat) : KM (List (Nat × AngleSum)) :=
  match angleSum #[beadLen beads s], mapK (nodeEntry beads Jt tree s) tree.toList with
  | .ok root, .ok rest => .ok ((s, root) :: rest)
  | .error e, _ => .error e
  | _, .error e => .error e

/-- The inequalities at one tree node `xs = (x, sx)` of the tree of `s`. -/
def checkNode (beads : Array Bead) (Jt : List (List (Nat × Bool)))
    (sums : List (Nat × AngleSum)) (s : Nat) (xs : Nat × AngleSum) : KM Unit :=
  forallK ((Jt.getD xs.1 []).filter (·.2)) fun yb =>
    if yb.1 = s then
      if atLeast xs.2 2 then .ok ()
      else .error s!"bead {s} closes a loop shorter than 2*pi through {xs.1}"
    else if yb.1 < s then .ok ()
    else
      match angleAdd xs.2 (beadLen beads yb.1) with
      | .error e => .error e
      | .ok sxy =>
        if atLeast sxy 2 then .ok ()
        else match sums.lookup yb.1 with
          | none => .error s!"tree of {s}: {yb.1} reachable via {xs.1} but missing"
          | some sy =>
            match angleLeq sy sxy with
            | .ok true => .ok ()
            | .ok false => .error s!"tree of {s}: distance to {yb.1} not a lower bound via {xs.1}"
            | .error e => .error e

/-- Verify the shortest-path certificate for one vertex link. -/
def checkVertexTrees (L : Link) (beads : Array Bead)
    (trees : Std.HashMap Nat (Std.HashMap Nat Nat)) : KM Unit :=
  match mapK (junctionRow L beads) (List.range beads.size) with
  | .error e => .error e
  | .ok Jt =>
    forallK (List.range beads.size) fun s =>
      match nodeSums beads Jt (trees.getD s {}) s with
      | .error e => .error e
      | .ok sums => forallK sums (checkNode beads Jt sums s)

lemma lookup_mem {l : List (Nat × AngleSum)} {x : Nat} {st : AngleSum}
    (h : l.lookup x = some st) : (x, st) ∈ l := by
  induction l with
  | nil => simp at h
  | cons p l ih =>
    obtain ⟨y, sy⟩ := p
    by_cases hxy : x = y
    · subst hxy
      simp at h
      subst h
      exact List.mem_cons_self
    · have : (x == y) = false := by simp [hxy]
      simp only [List.lookup, this] at h
      exact List.mem_cons_of_mem _ (ih h)

/-- The nodes of the tree of `s` with their sums (empty if the tree failed). -/
def sumsOf (beads : Array Bead) (Jt : List (List (Nat × Bool)))
    (trees : Std.HashMap Nat (Std.HashMap Nat Nat)) (s : Nat) : List (Nat × AngleSum) :=
  match nodeSums beads Jt (trees.getD s {}) s with
  | .ok l => l
  | .error _ => []

lemma nodeSums_ok {beads : Array Bead} {Jt : List (List (Nat × Bool))}
    {tree : Std.HashMap Nat Nat} {s : Nat} {sums : List (Nat × AngleSum)}
    (h : nodeSums beads Jt tree s = .ok sums) :
    (∃ root, sums.lookup s = some root ∧ root.Rep (angleOf (beadLen beads s).1 (beadLen beads s).2)) ∧
      ∀ p ∈ sums, p.2.Rep p.2.value := by
  unfold nodeSums at h
  split at h
  · rename_i root rest hroot hrest
    simp only [Except.ok.injEq] at h
    subst h
    have h2 := angleSum_ok hroot
    simp at h2
    refine ⟨⟨root, by simp, h2⟩, ?_⟩
    intro p hp
    rcases List.mem_cons.mp hp with rfl | hp
    · rw [AngleSum.value_eq h2]; exact h2
    · obtain ⟨a, -, ha⟩ := forall2_mem_right (mapK_ok hrest) hp
      obtain ⟨Θ, hΘ⟩ := nodeEntry_rep ha
      rw [AngleSum.value_eq hΘ]; exact hΘ
  · cases h
  · cases h

/-- Soundness of True certificates at one vertex: when the tree check
passes, every necklace — a cyclic sequence of beads, each passing to the
next — has total length at least `2π`. -/
theorem checkVertexTrees_ok {L : Link} {beads : Array Bead}
    {trees : Std.HashMap Nat (Std.HashMap Nat Nat)} (h : checkVertexTrees L beads trees = .ok ()) :
    ∃ Jt : List (List (Nat × Bool)),
      (∀ x < beads.size, junctionRow L beads x = .ok (Jt.getD x [])) ∧
      ∀ l : List ℕ, l ≠ [] → (∀ i ∈ l, i < beads.size) →
        (l : Cycle ℕ).Chain (fun x y => (y, true) ∈ Jt.getD x []) →
        2 * Real.pi ≤ (l.map fun i => angleOf (beadLen beads i).1 (beadLen beads i).2).sum := by
  unfold checkVertexTrees at h
  split at h
  · cases h
  rename_i Jt hJt
  refine ⟨Jt, ?_, ?_⟩
  · intro x hx
    have hf := mapK_ok hJt
    have hlen := hf.length_eq
    simp only [List.length_range] at hlen
    rw [List.getD_eq_getElem _ _ (by omega)]
    have := (List.forall₂_iff_get.mp hf).2 x (by simp [hx]) (by omega)
    simpa using this
  intro l hl hn hc
  have hall := forallK_ok h
  have hs : ∀ s < beads.size, nodeSums beads Jt (trees.getD s {}) s = .ok (sumsOf beads Jt trees s) ∧
      forallK (sumsOf beads Jt trees s) (checkNode beads Jt (sumsOf beads Jt trees s) s) = .ok () := by
    intro s hsn
    have := hall s (List.mem_range.mpr hsn)
    split at this
    · cases this
    · rename_i l' hl'
      simp only [sumsOf, hl']
      exact ⟨trivial, this⟩
  have hd : ∀ s < beads.size, ∀ x st, (sumsOf beads Jt trees s).lookup x = some st →
      st.Rep (((sumsOf beads Jt trees s).lookup x).elim 0 AngleSum.value) ∧
        (x, st) ∈ sumsOf beads Jt trees s := by
    intro s hsn x st hx
    have hm := lookup_mem hx
    refine ⟨?_, hm⟩
    simp only [hx, Option.elim]
    exact (nodeSums_ok (hs s hsn).1).2 _ hm
  refine no_short_necklace (n := beads.size)
    (len := fun i => angleOf (beadLen beads i).1 (beadLen beads i).2)
    (pass := fun x y => (y, true) ∈ Jt.getD x [])
    (T := fun s x => ((sumsOf beads Jt trees s).lookup x).isSome)
    (d := fun s x => ((sumsOf beads Jt trees s).lookup x).elim 0 AngleSum.value)
    (fun i => angleOf_nonneg _ _) ?_ ?_ ?_ l hl hn hc
  · -- the root
    intro s hsn
    obtain ⟨root, hlk, hR⟩ := (nodeSums_ok (hs s hsn).1).1
    refine ⟨by simp [hlk], ?_⟩
    rw [AngleSum.rep_unique (hd s hsn s root hlk).1 hR]
  · -- a step to a bead y > s
    intro s hsn x hT y hpass hsy hyn
    obtain ⟨st, hst⟩ := Option.isSome_iff_exists.mp hT
    obtain ⟨hR, hmem⟩ := hd s hsn x st hst
    have hcn := forallK_ok (hs s hsn).2 (x, st) hmem
    have hcheck := forallK_ok hcn (y, true) (List.mem_filter.mpr ⟨hpass, rfl⟩)
    simp only at hcheck
    rw [if_neg (by omega), if_neg (by omega)] at hcheck
    split at hcheck
    · cases hcheck
    · rename_i sxy hsxy
      have hR' := AngleSum.add_sound hR (angleAdd_ok hsxy)
      split_ifs at hcheck with hat
      · left
        have := (AngleSum.atLeast_iff hR' 2).mp (by simpa [atLeast] using hat)
        push_cast at this
        exact this
      · right
        split at hcheck
        · cases hcheck
        · rename_i sy hsy
          refine ⟨by simp [hsy], ?_⟩
          split at hcheck
          · rename_i hle
            exact (AngleSum.leq_sound (hd s hsn y sy hsy).1 hR' (angleLeq_ok hle)).mp rfl
          · cases hcheck
          · cases hcheck
  · -- closing the loop at s
    intro s hsn x hT hpass
    obtain ⟨st, hst⟩ := Option.isSome_iff_exists.mp hT
    obtain ⟨hR, hmem⟩ := hd s hsn x st hst
    have hcn := forallK_ok (hs s hsn).2 (x, st) hmem
    have hcheck := forallK_ok hcn (s, true) (List.mem_filter.mpr ⟨hpass, rfl⟩)
    simp only at hcheck
    rw [if_pos trivial] at hcheck
    split_ifs at hcheck with hat
    have := (AngleSum.atLeast_iff hR 2).mp (by simpa [atLeast] using hat)
    push_cast at this
    exact this

/-- A certificate that the metric is locally CAT(0) (`vlink-certificate-true/2`). -/
structure TrueCert where
  isosig : String
  sqlens : Array ℚ
  piBeads : Bool
  allowRepeats : Bool
  assumptions : Array String
  /-- per vertex: number of beads and the trees -/
  vertices : Array (Nat × Nat × Std.HashMap Nat (Std.HashMap Nat Nat))

/-- Every edge has dihedral sum `≥ 2π`. -/
def edgesOkK (T : Tri) (geo : Geometry) : KM Unit :=
  forallK (List.range T.edges.size) fun e =>
    match geo.edgeAngle e with
    | .ok acc => if atLeast acc 2 then .ok () else .error s!"edge {e} has dihedral sum < 2*pi"
    | .error err => .error err

/-- One vertex of a True certificate: the kernel's own bead enumeration, and
the shortest-path trees. -/
def vertexOkK (geo : Geometry) (piBeads : Bool)
    (v : Nat × Nat × Std.HashMap Nat (Std.HashMap Nat Nat)) : KM Unit :=
  match Link.make geo v.1 with
  | .error e => .error e
  | .ok L =>
    match enumerateBeads L piBeads with
    | .error e => .error e
    | .ok beads =>
      if v.2.1 ≠ beads.size then
        .error s!"vertex {v.1}: bead count differs from the kernel enumeration"
      else checkVertexTrees L beads v.2.2

def assumptionsNote (a : Array String) : String :=
  if a.isEmpty then "" else s!" (assuming {"; ".intercalate a.toList})"

/-- Check a certificate that the metric is locally CAT(0). -/
def checkTrue (c : TrueCert) : KM String :=
  if !c.piBeads && !c.assumptions.contains LEMMA_4_9 then
    .error s!"pi_beads=False requires the assumption: {LEMMA_4_9}"
  else if !c.allowRepeats then .error "allow_repeats=False relies on an unproven rule"
  else match Tri.ofIsoSig c.isosig with
  | .error e => .error e
  | .ok T =>
  match Geometry.make T c.sqlens with
  | .error e => .error e
  | .ok geo =>
  match tetsValidK T geo (fun t => s!"tetrahedron {t} is degenerate") with
  | .error e => .error e
  | .ok () =>
  match edgesOkK T geo with
  | .error e => .error e
  | .ok () =>
    let vs := c.vertices.qsort (fun a b => a.1 < b.1)
    if vs.toList.map (·.1) ≠ List.range T.corners.size then .error "certificate must cover every vertex"
    else match forallK vs.toList (vertexOkK geo c.piBeads) with
      | .error e => .error e
      | .ok () =>
        .ok s!"locally CAT(0): {T.n} tetrahedra valid, all edges >= 2*pi, no short necklace at {vs.size} vertices{assumptionsNote c.assumptions}"

end VLink
