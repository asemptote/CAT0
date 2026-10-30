import VLink.Angle

/-!
# Vectors of a developed vertex link (layer L3): exact predicates and their meaning

A developed strip of link triangles lives in `ℝ³` with the inner product
`⟨u, v⟩ = Σ G i j u_i v_j` of its first triangle's Gram matrix `G`; points of
the link are rays. `V3` holds a vector with exact multiquadratic coordinates.

For each predicate the kernel uses, this file proves what the exact
computation decides, as a statement about the real vectors:

* `val_dot`, `val_det3`, `val_cross`: the exact operations compute the inner
  product, determinant and cross product;
* `parallel_ok`: `parallel` decides whether two vectors are parallel (cross
  product zero) and, if so, the sign of their inner product;
* `crosses_ok`: when `crosses` accepts, the geodesic from `S` toward `Q`
  meets the open edge arc `AB` in its interior: some `αA + βB` (`α, β > 0`)
  equals `γS + δQ` with `γ, δ > 0` (minor arc), or with `δ > 0` (half great
  circle, for a symmetric `G`);
* `length_angle`, `tangent_angle`: the angle the kernel adds for a bead
  length is the angle between the two vectors, and for a corner the angle
  between the tangent directions at `W` (the projections orthogonal to `W`).
-/

namespace VLink

abbrev KM := Except String

def kfuel : Nat := 64

def sgn (x : MQ) : KM Int :=
  match MQ.sign kfuel x with
  | some s => .ok s
  | none => .error "sign undecided (out of fuel)"

def isZero (x : MQ) : KM Bool :=
  match MQ.sign kfuel x with
  | some s => .ok (s == 0)
  | none => .error "sign undecided (out of fuel)"

lemma sgn_ok {x : MQ} {s : Int} (h : sgn x = .ok s) : s = rsign x.val := by
  unfold sgn at h
  split at h
  · rename_i t ht
    cases h
    exact MQ.sign_sound _ _ _ ht
  · cases h

lemma isZero_ok {x : MQ} {b : Bool} (h : isZero x = .ok b) : (b = true ↔ x.val = 0) := by
  unfold isZero at h
  split at h
  · rename_i t ht
    cases h
    rw [MQ.sign_sound _ _ _ ht, beq_iff_eq, rsign_eq_zero]
  · cases h

/-- A vector with exact coordinates. -/
structure V3 where
  x : MQ
  y : MQ
  z : MQ
  deriving Inhabited

abbrev Gram := Array (Array ℚ)

namespace V3

def get (u : V3) : Nat → MQ
  | 0 => u.x
  | 1 => u.y
  | _ => u.z

def add (u v : V3) : V3 := ⟨MQ.add u.x v.x, MQ.add u.y v.y, MQ.add u.z v.z⟩
def sub (u v : V3) : V3 := ⟨MQ.sub u.x v.x, MQ.sub u.y v.y, MQ.sub u.z v.z⟩
def neg (u : V3) : V3 := ⟨MQ.neg u.x, MQ.neg u.y, MQ.neg u.z⟩
def scale (c : MQ) (u : V3) : V3 := ⟨MQ.mul c u.x, MQ.mul c u.y, MQ.mul c u.z⟩
def scaleQ (c : ℚ) (u : V3) : V3 := ⟨MQ.smul c u.x, MQ.smul c u.y, MQ.smul c u.z⟩

def cross (u v : V3) : V3 :=
  ⟨MQ.sub (MQ.mul u.y v.z) (MQ.mul u.z v.y), MQ.sub (MQ.mul u.z v.x) (MQ.mul u.x v.z),
   MQ.sub (MQ.mul u.x v.y) (MQ.mul u.y v.x)⟩

noncomputable def val (u : V3) : Fin 3 → ℝ := ![u.x.val, u.y.val, u.z.val]

lemma val_apply (u : V3) :
    u.val 0 = u.x.val ∧ u.val 1 = u.y.val ∧ u.val 2 = u.z.val := ⟨rfl, rfl, rfl⟩

@[simp] lemma val_add (u v : V3) : (add u v).val = u.val + v.val := by
  ext i; fin_cases i <;> simp [add, val, MQ.val_add]

@[simp] lemma val_sub (u v : V3) : (sub u v).val = u.val - v.val := by
  ext i; fin_cases i <;> simp [sub, val, MQ.val_sub]

@[simp] lemma val_neg (u : V3) : (neg u).val = -u.val := by
  ext i; fin_cases i <;> simp [neg, val, MQ.val_neg]

@[simp] lemma val_scale (c : MQ) (u : V3) : (scale c u).val = c.val • u.val := by
  ext i; fin_cases i <;> simp [scale, val, MQ.val_mul]

@[simp] lemma val_scaleQ (c : ℚ) (u : V3) : (scaleQ c u).val = (c : ℝ) • u.val := by
  ext i; fin_cases i <;> simp [scaleQ, val, MQ.val_smul]

lemma val_cross (u v : V3) : (cross u v).val = crossProduct u.val v.val := by
  rw [cross_apply]
  ext i; fin_cases i <;> simp [cross, val, MQ.val_sub, MQ.val_mul]

end V3

open V3

/-! ### Inner product and determinant -/

/-- The Gram matrix as a real matrix. -/
noncomputable def Gram.toMat (G : Gram) : Fin 3 → Fin 3 → ℝ :=
  fun i j => (((G[i.val]!)[j.val]! : ℚ) : ℝ)

/-- `⟨u, v⟩ = Σ M i j u_i v_j`. -/
noncomputable def bform (M : Fin 3 → Fin 3 → ℝ) (u v : Fin 3 → ℝ) : ℝ :=
  ∑ i, ∑ j, M i j * u i * v j

def dotTerms : List (Nat × Nat) :=
  [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]

def dotStep (G : Gram) (u v : V3) (acc : MQ) (jk : Nat × Nat) : MQ :=
  let g := (G[jk.1]!)[jk.2]!
  if g = 0 then acc else MQ.add acc (MQ.smul g (MQ.mul (u.get jk.1) (v.get jk.2)))

def dot (G : Gram) (u v : V3) : MQ := dotTerms.foldl (dotStep G u v) []

lemma dot_fold (G : Gram) (u v : V3) (l : List (Nat × Nat)) (acc : MQ) :
    (l.foldl (dotStep G u v) acc).val =
      acc.val + (l.map fun jk => (((G[jk.1]!)[jk.2]! : ℚ) : ℝ) *
        (u.get jk.1).val * (v.get jk.2).val).sum := by
  induction l generalizing acc with
  | nil => simp
  | cons jk l ih =>
    rw [List.foldl_cons, ih, List.map_cons, List.sum_cons]
    unfold dotStep
    dsimp only
    split_ifs with h
    · rw [h]; simp
    · rw [MQ.val_add, MQ.val_smul, MQ.val_mul]; ring

lemma val_dot (G : Gram) (u v : V3) : (dot G u v).val = bform G.toMat u.val v.val := by
  rw [dot, dot_fold]
  simp [dotTerms, bform, Fin.sum_univ_three, Gram.toMat, V3.val, V3.get]
  ring

def det3 (u v w : V3) : MQ :=
  let m := MQ.mul
  let s := MQ.sub
  MQ.add (s (m u.x (s (m v.y w.z) (m v.z w.y))) (m u.y (s (m v.x w.z) (m v.z w.x))))
    (m u.z (s (m v.x w.y) (m v.y w.x)))

/-- The real determinant of the rows `u, v, w`. -/
noncomputable def rdet (u v w : Fin 3 → ℝ) : ℝ := Matrix.det (Matrix.of ![u, v, w])

lemma rdet_eq (u v w : Fin 3 → ℝ) :
    rdet u v w = u 0 * (v 1 * w 2 - v 2 * w 1) - u 1 * (v 0 * w 2 - v 2 * w 0) +
      u 2 * (v 0 * w 1 - v 1 * w 0) := by
  rw [rdet, Matrix.det_fin_three]
  simp only [Matrix.of_apply, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
    Matrix.head_cons, Matrix.tail_cons]
  ring

lemma val_det3 (u v w : V3) : (det3 u v w).val = rdet u.val v.val w.val := by
  rw [rdet_eq]
  simp [det3, MQ.val_add, MQ.val_sub, MQ.val_mul, V3.val]

/-- Cramer: `det(A,B,C) X = det(X,B,C) A + det(A,X,C) B + det(A,B,X) C`. -/
lemma cramer (A B C X : Fin 3 → ℝ) :
    rdet A B C • X = rdet X B C • A + rdet A X C • B + rdet A B X • C := by
  ext i
  simp only [rdet_eq, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  fin_cases i <;> simp <;> ring

/-! ### Parallel vectors -/

/-- `+1` same direction, `-1` opposite, `0` not parallel. -/
def parallel (G : Gram) (x y : V3) : KM Int :=
  match isZero (V3.cross x y).x, isZero (V3.cross x y).y, isZero (V3.cross x y).z with
  | .ok zx, .ok zy, .ok zz => if zx && zy && zz then sgn (dot G x y) else .ok 0
  | _, _, _ => .error "sign undecided (out of fuel)"

theorem parallel_ok {G : Gram} {x y : V3} {r : Int} (h : parallel G x y = .ok r) :
    (crossProduct x.val y.val = 0 ∧ r = rsign (bform G.toMat x.val y.val)) ∨
      (crossProduct x.val y.val ≠ 0 ∧ r = 0) := by
  unfold parallel at h
  split at h
  · rename_i zx zy zz hx hy hz
    have ex := isZero_ok hx
    have ey := isZero_ok hy
    have ez := isZero_ok hz
    have hc : crossProduct x.val y.val =
        ![(V3.cross x y).x.val, (V3.cross x y).y.val, (V3.cross x y).z.val] := by
      rw [← val_cross]; rfl
    split_ifs at h with hall
    · left
      simp only [Bool.and_eq_true] at hall
      refine ⟨?_, ?_⟩
      · rw [hc, ex.mp hall.1.1, ey.mp hall.1.2, ez.mp hall.2]
        ext i; fin_cases i <;> simp
      · rw [sgn_ok h, val_dot]
    · right
      refine ⟨?_, by cases h; rfl⟩
      intro h0
      rw [hc] at h0
      have h0x := congrFun h0 0
      have h0y := congrFun h0 1
      have h0z := congrFun h0 2
      simp at h0x h0y h0z
      apply hall
      simp [ex.mpr h0x, ey.mpr h0y, ez.mpr h0z]
  · cases h

/-! ### Geodesics crossing edges -/

def coneSigns (X A B C : V3) : KM (Int × Int) :=
  match sgn (det3 A B C), sgn (det3 X B C), sgn (det3 A X C) with
  | .ok den, .ok a, .ok b => .ok (a * den, b * den)
  | _, _, _ => .error "sign undecided (out of fuel)"

lemma coneSigns_ok {X A B C : V3} {al be : Int} (h : coneSigns X A B C = .ok (al, be)) :
    al = rsign (rdet X.val B.val C.val) * rsign (rdet A.val B.val C.val) ∧
      be = rsign (rdet A.val X.val C.val) * rsign (rdet A.val B.val C.val) := by
  unfold coneSigns at h
  split at h
  · rename_i den a b hd ha hb
    simp only [Except.ok.injEq, Prod.mk.injEq] at h
    obtain ⟨rfl, rfl⟩ := h
    rw [sgn_ok ha, sgn_ok hb, sgn_ok hd, val_det3, val_det3, val_det3]
    exact ⟨rfl, rfl⟩
  · cases h

/-- Does the geodesic from `S` toward `Q` cross the open edge arc `AB`?
`halfCircle = false`: the minor arc from `S` to `Q`; `true`: the half great
circle from `S` through `Q`'s side to `-S`. -/
def crosses (G : Gram) (S Q A B : V3) (halfCircle : Bool) : KM Bool :=
  match isZero (det3 A B S) with
  | .error e => .error e
  | .ok true => .ok false
  | .ok false =>
    let X := V3.sub (V3.scale (det3 A B Q) S) (V3.scale (det3 A B S) Q)
    match isZero X.x, isZero X.y, isZero X.z with
    | .ok zx, .ok zy, .ok zz =>
      if zx && zy && zz then .ok false
      else
        match coneSigns X A B S with
        | .error e => .error e
        | .ok (al, be) =>
          if al = be ∧ al ≠ 0 then
            let Y := if al > 0 then X else V3.neg X
            if halfCircle then
              match sgn (MQ.sub (MQ.mul (dot G Y Q) (dot G S S))
                  (MQ.mul (dot G Y S) (dot G Q S))) with
              | .ok s => .ok (decide (s > 0))
              | .error e => .error e
            else
              match isZero (det3 S Q A) with
              | .error e => .error e
              | .ok zA =>
                match coneSigns Y S Q (if zA then B else A) with
                | .ok (ga, de) => .ok (decide (ga > 0 ∧ de > 0))
                | .error e => .error e
          else .ok false
    | _, _, _ => .error "sign undecided (out of fuel)"

lemma rsign_prod_pos {a b : ℝ} (h : 0 < rsign a * rsign b) : 0 < a * b := by
  rcases rsign_cases a with ⟨ha, ea⟩ | ⟨ha, ea⟩ | ⟨ha, ea⟩ <;>
    rcases rsign_cases b with ⟨hb, eb⟩ | ⟨hb, eb⟩ | ⟨hb, eb⟩ <;>
    rw [ea, eb] at h <;> first | nlinarith | (norm_num at h)

lemma rsign_prod_neg {a b : ℝ} (h : rsign a * rsign b < 0) : a * b < 0 := by
  rcases rsign_cases a with ⟨ha, ea⟩ | ⟨ha, ea⟩ | ⟨ha, ea⟩ <;>
    rcases rsign_cases b with ⟨hb, eb⟩ | ⟨hb, eb⟩ | ⟨hb, eb⟩ <;>
    rw [ea, eb] at h <;> first | nlinarith | (norm_num at h)

/-- A vector in the plane of `P₁, P₂` whose Cramer coefficients (against any
`C` off that plane) are positive is a positive combination of `P₁, P₂`. -/
lemma cone_of_signs {P1 P2 C Y : Fin 3 → ℝ} (h0 : rdet P1 P2 Y = 0)
    (h1 : 0 < rdet Y P2 C * rdet P1 P2 C) (h2 : 0 < rdet P1 Y C * rdet P1 P2 C) :
    ∃ γ δ : ℝ, 0 < γ ∧ 0 < δ ∧ Y = γ • P1 + δ • P2 := by
  have hD : rdet P1 P2 C ≠ 0 := by
    intro h; rw [h, mul_zero] at h1; exact lt_irrefl _ h1
  have hD2 : 0 < rdet P1 P2 C ^ 2 := by positivity
  refine ⟨rdet Y P2 C / rdet P1 P2 C, rdet P1 Y C / rdet P1 P2 C, ?_, ?_, ?_⟩
  · have e : rdet Y P2 C / rdet P1 P2 C = rdet Y P2 C * rdet P1 P2 C / rdet P1 P2 C ^ 2 := by
      field_simp
    rw [e]; exact div_pos h1 hD2
  · have e : rdet P1 Y C / rdet P1 P2 C = rdet P1 Y C * rdet P1 P2 C / rdet P1 P2 C ^ 2 := by
      field_simp
    rw [e]; exact div_pos h2 hD2
  · have cr := cramer P1 P2 C Y
    rw [h0, zero_smul, add_zero] at cr
    apply smul_right_injective (Fin 3 → ℝ) hD
    simp only
    rw [cr, smul_add, smul_smul, smul_smul]
    congr 2 <;> field_simp

lemma rdet_smul_left (c : ℝ) (X B C : Fin 3 → ℝ) : rdet (c • X) B C = c * rdet X B C := by
  simp only [rdet_eq, Pi.smul_apply, smul_eq_mul]; ring

lemma rdet_smul_mid (c : ℝ) (A X C : Fin 3 → ℝ) : rdet A (c • X) C = c * rdet A X C := by
  simp only [rdet_eq, Pi.smul_apply, smul_eq_mul]; ring

lemma rdet_smul_right (c : ℝ) (A B X : Fin 3 → ℝ) : rdet A B (c • X) = c * rdet A B X := by
  simp only [rdet_eq, Pi.smul_apply, smul_eq_mul]; ring

/-- `X = det(A,B,Q) S - det(A,B,S) Q` lies on both planes: of `A, B` and of `S, Q`. -/
lemma rdet_X (A B S Q : Fin 3 → ℝ) :
    rdet A B (rdet A B Q • S - rdet A B S • Q) = 0 ∧
      rdet S Q (rdet A B Q • S - rdet A B S • Q) = 0 := by
  constructor <;> simp only [rdet_eq, Pi.sub_apply, Pi.smul_apply, smul_eq_mul] <;> ring

lemma int_sign_cases (s : Int) (h : ∃ a b : ℝ, s = rsign a * rsign b) :
    s = 1 ∨ s = 0 ∨ s = -1 := by
  obtain ⟨a, b, rfl⟩ := h
  rcases rsign_cases a with ⟨_, ea⟩ | ⟨_, ea⟩ | ⟨_, ea⟩ <;>
    rcases rsign_cases b with ⟨_, eb⟩ | ⟨_, eb⟩ | ⟨_, eb⟩ <;> rw [ea, eb] <;> simp

/-- The first half of `crosses`: `Y = ±X` is a positive combination of `A, B`. -/
lemma crosses_cone_AB {A B S X : Fin 3 → ℝ} {al be : Int}
    (hal : al = rsign (rdet X B S) * rsign (rdet A B S))
    (hbe : be = rsign (rdet A X S) * rsign (rdet A B S))
    (hX : rdet A B X = 0) (hab : al = be) (hne : al ≠ 0) :
    ∃ α β : ℝ, 0 < α ∧ 0 < β ∧
      (if al > 0 then X else -X) = α • A + β • B := by
  by_cases hp : al > 0
  · rw [if_pos hp]
    refine cone_of_signs hX (rsign_prod_pos (hal ▸ hp)) (rsign_prod_pos ((hab ▸ hbe) ▸ hp))
  · rw [if_neg hp]
    have hn : al < 0 := lt_of_le_of_ne (not_lt.mp hp) hne
    have h1 := rsign_prod_neg (hal ▸ hn)
    have h2 := rsign_prod_neg ((hab ▸ hbe) ▸ hn)
    have e1 : rdet (-X) B S = -rdet X B S := by
      rw [show -X = (-1 : ℝ) • X by simp, rdet_smul_left]; ring
    have e2 : rdet A (-X) S = -rdet A X S := by
      rw [show -X = (-1 : ℝ) • X by simp, rdet_smul_mid]; ring
    have e0 : rdet A B (-X) = 0 := by
      rw [show -X = (-1 : ℝ) • X by simp, rdet_smul_right, hX, mul_zero]
    exact cone_of_signs e0 (by rw [e1]; linarith) (by rw [e2]; linarith)

/-- Soundness of `crosses` for the minor arc: the geodesic from `S` to `Q`
meets the open arc `AB` at an interior point of both. -/
theorem crosses_ok_minor {G : Gram} {S Q A B : V3}
    (h : crosses G S Q A B false = .ok true) :
    ∃ α β γ δ : ℝ, 0 < α ∧ 0 < β ∧ 0 < γ ∧ 0 < δ ∧
      α • A.val + β • B.val = γ • S.val + δ • Q.val := by
  unfold crosses at h
  simp only [Bool.false_eq_true, ↓reduceIte] at h
  cases hzS : isZero (det3 A B S) with
  | error e => simp [hzS] at h
  | ok zS =>
  cases zS with
  | true => simp [hzS] at h
  | false =>
  simp only [hzS] at h
  set X := V3.sub (V3.scale (det3 A B Q) S) (V3.scale (det3 A B S) Q) with hXdef
  cases hzx : isZero X.x with
  | error e => simp [hzx] at h
  | ok zx =>
  cases hzy : isZero X.y with
  | error e => simp [hzx, hzy] at h
  | ok zy =>
  cases hzz : isZero X.z with
  | error e => simp [hzx, hzy, hzz] at h
  | ok zz =>
  simp only [hzx, hzy, hzz] at h
  by_cases hall : (zx && zy && zz) = true
  · simp [hall] at h
  rw [if_neg hall] at h
  cases hcs : coneSigns X A B S with
  | error e => simp [hcs] at h
  | ok albe =>
  obtain ⟨al, be⟩ := albe
  simp only [hcs] at h
  by_cases hab : al = be ∧ al ≠ 0
  swap
  · rw [if_neg hab] at h; simp at h
  rw [if_pos hab] at h
  set Y := if al > 0 then X else V3.neg X with hYdef
  cases hzA : isZero (det3 S Q A) with
  | error e => simp [hzA] at h
  | ok zA =>
  simp only [hzA] at h
  cases hcs2 : coneSigns Y S Q (if zA = true then B else A) with
  | error e => simp [hcs2] at h
  | ok gade =>
  obtain ⟨ga, de⟩ := gade
  simp only [hcs2] at h
  have hpos : ga > 0 ∧ de > 0 := by simpa using h
  obtain ⟨hal, hbe⟩ := coneSigns_ok hcs
  have hXval : X.val = (det3 A B Q).val • S.val - (det3 A B S).val • Q.val := by
    rw [hXdef]; simp
  rw [val_det3, val_det3] at hXval
  obtain ⟨hXAB, hXSQ⟩ := rdet_X A.val B.val S.val Q.val
  rw [← hXval] at hXAB hXSQ
  obtain ⟨α, β, hα, hβ, hY⟩ := crosses_cone_AB hal hbe hXAB hab.1 hab.2
  have hYval : Y.val = if al > 0 then X.val else -X.val := by
    rw [hYdef]; split_ifs <;> simp
  rw [← hYval] at hY
  have hYSQ : rdet S.val Q.val Y.val = 0 := by
    rw [hYval]; split_ifs
    · exact hXSQ
    · rw [show -X.val = (-1 : ℝ) • X.val by simp, rdet_smul_right, hXSQ, mul_zero]
  obtain ⟨hga, hde⟩ := coneSigns_ok hcs2
  obtain ⟨γ, δ, hγ, hδ, hY2⟩ := cone_of_signs hYSQ
    (rsign_prod_pos (hga ▸ hpos.1)) (rsign_prod_pos (hde ▸ hpos.2))
  exact ⟨α, β, γ, δ, hα, hβ, hγ, hδ, hY.symm.trans hY2⟩

/-! ### Angles -/

/-- `(a, b²)` of the angle at `W` between the geodesics toward `X` and `Y`. -/
def tangentAngle (G : Gram) (W X Y : V3) : MQ × MQ :=
  let ww := dot G W W
  let xw := dot G X W
  let yw := dot G Y W
  let a := MQ.sub (MQ.mul (dot G X Y) ww) (MQ.mul xw yw)
  let nx := MQ.sub (MQ.mul (dot G X X) ww) (MQ.mul xw xw)
  let ny := MQ.sub (MQ.mul (dot G Y Y) ww) (MQ.mul yw yw)
  (a, MQ.sub (MQ.mul nx ny) (MQ.mul a a))

/-- `(a, b²)` of the angle between `S` and `E`. -/
def length (G : Gram) (S E : V3) : MQ × MQ :=
  let a := dot G S E
  (a, MQ.sub (MQ.mul (dot G S S) (dot G E E)) (MQ.mul a a))

lemma bform_sub_smul (M : Fin 3 → Fin 3 → ℝ) (X Y W : Fin 3 → ℝ) (c d : ℝ) :
    bform M (X - c • W) (Y - d • W) =
      bform M X Y - d * bform M X W - c * bform M W Y + c * d * bform M W W := by
  simp only [bform, Fin.sum_univ_three, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]; ring

lemma bform_add_smul (M : Fin 3 → Fin 3 → ℝ) (S Q Z : Fin 3 → ℝ) (c d : ℝ) :
    bform M (c • S + d • Q) Z = c * bform M S Z + d * bform M Q Z := by
  simp only [bform, Fin.sum_univ_three, Pi.add_apply, Pi.smul_apply, smul_eq_mul]; ring

/-- A bead's length: the angle the kernel adds for the bead from `S` to `E`
is the angle between `S` and `E`. -/
theorem length_angle {G : Gram} {S E : V3}
    (hb : 0 ≤ (length G S E).2.val) (hne : (length G S E).1.val ≠ 0 ∨ 0 < (length G S E).2.val)
    (hS : 0 ≤ bform G.toMat S.val S.val) :
    angleOf (length G S E).1 (length G S E).2 =
      Real.arccos (bform G.toMat S.val E.val /
        (Real.sqrt (bform G.toMat S.val S.val) * Real.sqrt (bform G.toMat E.val E.val))) := by
  rw [angleOf_eq_arccos _ _ hb hne]
  have h1 : (length G S E).1.val = bform G.toMat S.val E.val := by simp [length, val_dot]
  have h2 : (length G S E).2.val = bform G.toMat S.val S.val * bform G.toMat E.val E.val -
      bform G.toMat S.val E.val ^ 2 := by
    simp [length, MQ.val_sub, MQ.val_mul, val_dot]; ring
  rw [h1, h2, show ∀ a b : ℝ, a ^ 2 + (b - a ^ 2) = b by intros; ring, Real.sqrt_mul hS]

/-- A corner: the angle the kernel adds at `W` between the geodesics toward
`X` and `Y` is the angle between their tangent directions `X₁`, `Y₁`, the
components of `X` and `Y` orthogonal to `W` (for a symmetric form). -/
theorem tangent_angle {G : Gram} {W X Y : V3}
    (hsym : ∀ u v, bform G.toMat u v = bform G.toMat v u)
    (hW : 0 < bform G.toMat W.val W.val)
    (hb : 0 ≤ (tangentAngle G W X Y).2.val)
    (hne : (tangentAngle G W X Y).1.val ≠ 0 ∨ 0 < (tangentAngle G W X Y).2.val)
    (X₁ Y₁ : Fin 3 → ℝ)
    (hX₁ : X₁ = X.val - (bform G.toMat X.val W.val / bform G.toMat W.val W.val) • W.val)
    (hY₁ : Y₁ = Y.val - (bform G.toMat Y.val W.val / bform G.toMat W.val W.val) • W.val)
    (hX : 0 ≤ bform G.toMat X₁ X₁) :
    angleOf (tangentAngle G W X Y).1 (tangentAngle G W X Y).2 =
      Real.arccos (bform G.toMat X₁ Y₁ /
        (Real.sqrt (bform G.toMat X₁ X₁) * Real.sqrt (bform G.toMat Y₁ Y₁))) := by
  rw [angleOf_eq_arccos _ _ hb hne]
  set M := G.toMat with hM
  set ww := bform M W.val W.val with hww
  set xw := bform M X.val W.val with hxw
  set yw := bform M Y.val W.val with hyw
  have hWY : bform M W.val Y.val = yw := hsym _ _
  have hWX : bform M W.val X.val = xw := hsym _ _
  have hw0 : ww ≠ 0 := hW.ne'
  have eXY : bform M X₁ Y₁ = (bform M X.val Y.val * ww - xw * yw) / ww := by
    rw [hX₁, hY₁, bform_sub_smul, hWY]; field_simp; ring
  have eXX : bform M X₁ X₁ = (bform M X.val X.val * ww - xw * xw) / ww := by
    rw [hX₁, bform_sub_smul, hWX]; field_simp; ring
  have eYY : bform M Y₁ Y₁ = (bform M Y.val Y.val * ww - yw * yw) / ww := by
    rw [hY₁, bform_sub_smul, hWY]; field_simp; ring
  have ha : (tangentAngle G W X Y).1.val = ww * bform M X₁ Y₁ := by
    rw [eXY]
    simp only [tangentAngle, MQ.val_sub, MQ.val_mul, val_dot]
    rw [← hM, ← hww, ← hxw, ← hyw]
    field_simp
  have hb2 : (tangentAngle G W X Y).2.val =
      (ww * bform M X₁ X₁) * (ww * bform M Y₁ Y₁) - (ww * bform M X₁ Y₁) ^ 2 := by
    rw [eXX, eYY, eXY]
    simp only [tangentAngle, MQ.val_sub, MQ.val_mul, val_dot]
    rw [← hM, ← hww, ← hxw, ← hyw]
    field_simp
  have hY : 0 ≤ bform M Y₁ Y₁ := by
    by_contra hneg
    push Not at hneg
    have hsum : (tangentAngle G W X Y).2.val + (tangentAngle G W X Y).1.val ^ 2 =
        (ww * bform M X₁ X₁) * (ww * bform M Y₁ Y₁) := by rw [hb2, ha]; ring
    rcases hX.lt_or_eq with hx | hx
    · have hp : 0 < ww * bform M X₁ X₁ * ww := mul_pos (mul_pos hW hx) hW
      nlinarith [sq_nonneg ((tangentAngle G W X Y).1.val), mul_neg_of_pos_of_neg hp hneg]
    · rw [← hx] at hsum
      rcases hne with hne | hne
      · nlinarith [sq_pos_of_ne_zero hne]
      · nlinarith
  rw [ha, hb2, show ∀ a b : ℝ, a ^ 2 + (b - a ^ 2) = b by intros; ring,
    Real.sqrt_mul (by positivity), Real.sqrt_mul hW.le, Real.sqrt_mul hW.le]
  rw [show Real.sqrt ww * Real.sqrt (bform M X₁ X₁) * (Real.sqrt ww * Real.sqrt (bform M Y₁ Y₁)) =
      ww * (Real.sqrt (bform M X₁ X₁) * Real.sqrt (bform M Y₁ Y₁)) by
    rw [show Real.sqrt ww * Real.sqrt (bform M X₁ X₁) * (Real.sqrt ww * Real.sqrt (bform M Y₁ Y₁)) =
      (Real.sqrt ww * Real.sqrt ww) * (Real.sqrt (bform M X₁ X₁) * Real.sqrt (bform M Y₁ Y₁)) by ring,
      Real.mul_self_sqrt hW.le]]
  rw [mul_div_mul_left _ _ hw0]

/-- Soundness of `crosses` for the half great circle (π-beads), for a
symmetric form: the half circle from `S` through the side of `Q`, the rays
`γS + δQ` with `δ > 0`, meets the open arc `AB`. -/
theorem crosses_ok_half {G : Gram} {S Q A B : V3}
    (hsym : ∀ u v, bform G.toMat u v = bform G.toMat v u)
    (h : crosses G S Q A B true = .ok true) :
    ∃ α β γ δ : ℝ, 0 < α ∧ 0 < β ∧
      α • A.val + β • B.val = γ • S.val + δ • Q.val ∧
      0 < δ * (bform G.toMat Q.val Q.val * bform G.toMat S.val S.val -
        bform G.toMat Q.val S.val ^ 2) := by
  unfold crosses at h
  simp only [↓reduceIte] at h
  cases hzS : isZero (det3 A B S) with
  | error e => simp [hzS] at h
  | ok zS =>
  cases zS with
  | true => simp [hzS] at h
  | false =>
  simp only [hzS] at h
  set X := V3.sub (V3.scale (det3 A B Q) S) (V3.scale (det3 A B S) Q) with hXdef
  cases hzx : isZero X.x with
  | error e => simp [hzx] at h
  | ok zx =>
  cases hzy : isZero X.y with
  | error e => simp [hzx, hzy] at h
  | ok zy =>
  cases hzz : isZero X.z with
  | error e => simp [hzx, hzy, hzz] at h
  | ok zz =>
  simp only [hzx, hzy, hzz] at h
  by_cases hall : (zx && zy && zz) = true
  · simp [hall] at h
  rw [if_neg hall] at h
  cases hcs : coneSigns X A B S with
  | error e => simp [hcs] at h
  | ok albe =>
  obtain ⟨al, be⟩ := albe
  simp only [hcs] at h
  by_cases hab : al = be ∧ al ≠ 0
  swap
  · rw [if_neg hab] at h; simp at h
  rw [if_pos hab] at h
  set Y := if al > 0 then X else V3.neg X with hYdef
  set e := MQ.sub (MQ.mul (dot G Y Q) (dot G S S)) (MQ.mul (dot G Y S) (dot G Q S)) with hedef
  cases hs : sgn e with
  | error err => simp [hs] at h
  | ok sv =>
  simp only [hs] at h
  have hsv : sv > 0 := by simpa using h
  have hepos : 0 < e.val := by
    rw [sgn_ok hs] at hsv
    rcases rsign_cases e.val with ⟨h', _⟩ | ⟨_, e'⟩ | ⟨_, e'⟩
    · exact h'
    · rw [e'] at hsv; norm_num at hsv
    · rw [e'] at hsv; norm_num at hsv
  obtain ⟨hal, hbe⟩ := coneSigns_ok hcs
  have hXval : X.val = (det3 A B Q).val • S.val - (det3 A B S).val • Q.val := by
    rw [hXdef]; simp
  rw [val_det3, val_det3] at hXval
  obtain ⟨hXAB, -⟩ := rdet_X A.val B.val S.val Q.val
  rw [← hXval] at hXAB
  obtain ⟨α, β, hα, hβ, hY⟩ := crosses_cone_AB hal hbe hXAB hab.1 hab.2
  have hYval : Y.val = if al > 0 then X.val else -X.val := by
    rw [hYdef]; split_ifs <;> simp
  rw [← hYval] at hY
  set σ : ℝ := if al > 0 then 1 else -1 with hσ
  set dQ := rdet A.val B.val Q.val
  set dS := rdet A.val B.val S.val
  have hYcomb : Y.val = (σ * dQ) • S.val + (-(σ * dS)) • Q.val := by
    rw [hYval, hXval, hσ]
    split_ifs
    · ext i; simp [Pi.sub_apply, Pi.smul_apply, smul_eq_mul]; ring
    · ext i; simp [Pi.sub_apply, Pi.smul_apply, smul_eq_mul]; ring
  refine ⟨α, β, σ * dQ, -(σ * dS), hα, hβ, hY.symm.trans hYcomb, ?_⟩
  have hev : e.val = bform G.toMat Y.val Q.val * bform G.toMat S.val S.val -
      bform G.toMat Y.val S.val * bform G.toMat Q.val S.val := by
    simp [hedef, MQ.val_sub, MQ.val_mul, val_dot]
  rw [hev, hYcomb, bform_add_smul, bform_add_smul, hsym S.val Q.val] at hepos
  nlinarith [hepos]

end VLink
