import VLink.Vec

/-!
# Developing a strip is isometric (layer L3)

Crossing from link triangle `t` to its neighbour `t'` across the edge
opposite vertex `d` (`Kernel.step`): with `Pa, Pb` the shared edge's vertices
and `Pd` the third vertex of `t`, the kernel writes `Pd = p₀ Pa + p₁ Pb + n`
with `n ⊥ Pa, Pb` (`p` from the Gram matrix of `t`), and places the new
vertex at `q₀ Pa + q₁ Pb - μ n` (`q` from the Gram matrix of `t'`,
`μ = √(h'/h)`, `h = |n|²`). `step_iso` proves that the new vertex has exactly
the inner products of `t'` with `Pa`, `Pb` and itself, and lies on the other
side of the plane of `Pa, Pb` from `Pd`: the strip is an isometric copy of
the triangles glued along their shared edges.
-/

namespace VLink

/-- Bilinearity: `⟨aX + bY + cZ, W⟩`. -/
lemma bform_comb_left (M : Fin 3 → Fin 3 → ℝ) (X Y Z W : Fin 3 → ℝ) (a b c : ℝ) :
    bform M (a • X + b • Y + c • Z) W = a * bform M X W + b * bform M Y W + c * bform M Z W := by
  simp only [bform, Fin.sum_univ_three, Pi.add_apply, Pi.smul_apply, smul_eq_mul]; ring

/-- Bilinearity: `⟨aX + bY + cZ, a'X + b'Y + c'Z⟩`. -/
lemma bform_comb (M : Fin 3 → Fin 3 → ℝ) (X Y Z : Fin 3 → ℝ) (a b c a' b' c' : ℝ) :
    bform M (a • X + b • Y + c • Z) (a' • X + b' • Y + c' • Z) =
      a * a' * bform M X X + a * b' * bform M X Y + a * c' * bform M X Z +
      b * a' * bform M Y X + b * b' * bform M Y Y + b * c' * bform M Y Z +
      c * a' * bform M Z X + c * b' * bform M Z Y + c * c' * bform M Z Z := by
  simp only [bform, Fin.sum_univ_three, Pi.add_apply, Pi.smul_apply, smul_eq_mul]; ring

/-- One step of the development. `gaa … gdd`: the Gram matrix of `t` on its
vertices `a, b, d`; `ea, eb, ed`: the inner products of `t'`'s new vertex with
the shared edge's ends and itself. -/
theorem step_iso (M : Fin 3 → Fin 3 → ℝ) (hsym : ∀ u v, bform M u v = bform M v u)
    (Pa Pb Pd : Fin 3 → ℝ) (gaa gab gbb gda gdb gdd ea eb ed p0 p1 q0 q1 μ : ℝ)
    (haa : bform M Pa Pa = gaa) (hab : bform M Pa Pb = gab) (hbb : bform M Pb Pb = gbb)
    (hda : bform M Pd Pa = gda) (hdb : bform M Pd Pb = gdb) (hdd : bform M Pd Pd = gdd)
    (hdt : gaa * gbb - gab ^ 2 ≠ 0)
    (hp0 : p0 = (gda * gbb - gdb * gab) / (gaa * gbb - gab ^ 2))
    (hp1 : p1 = (gdb * gaa - gda * gab) / (gaa * gbb - gab ^ 2))
    (hq0 : q0 = (ea * gbb - eb * gab) / (gaa * gbb - gab ^ 2))
    (hq1 : q1 = (eb * gaa - ea * gab) / (gaa * gbb - gab ^ 2))
    (hμ : μ ^ 2 * (gdd - p0 * gda - p1 * gdb) = ed - q0 * ea - q1 * eb) :
    let n := Pd - p0 • Pa - p1 • Pb
    let new := q0 • Pa + q1 • Pb - μ • n
    bform M new Pa = ea ∧ bform M new Pb = eb ∧ bform M new new = ed ∧
      rdet Pa Pb new = -μ * rdet Pa Pb Pd := by
  intro n new
  set D := gaa * gbb - gab ^ 2 with hD
  have ep0 : p0 * D = gda * gbb - gdb * gab := by rw [hp0, div_mul_cancel₀ _ hdt]
  have ep1 : p1 * D = gdb * gaa - gda * gab := by rw [hp1, div_mul_cancel₀ _ hdt]
  have eq0 : q0 * D = ea * gbb - eb * gab := by rw [hq0, div_mul_cancel₀ _ hdt]
  have eq1 : q1 * D = eb * gaa - ea * gab := by rw [hq1, div_mul_cancel₀ _ hdt]
  have I1 : gda - p0 * gaa - p1 * gab = 0 := by
    have : (gda - p0 * gaa - p1 * gab) * D = 0 := by
      linear_combination (-gaa) * ep0 + (-gab) * ep1 + gda * hD
    exact (mul_eq_zero.mp this).resolve_right hdt
  have I2 : gdb - p0 * gab - p1 * gbb = 0 := by
    have : (gdb - p0 * gab - p1 * gbb) * D = 0 := by
      linear_combination (-gab) * ep0 + (-gbb) * ep1 + gdb * hD
    exact (mul_eq_zero.mp this).resolve_right hdt
  have I3 : ea - q0 * gaa - q1 * gab = 0 := by
    have : (ea - q0 * gaa - q1 * gab) * D = 0 := by
      linear_combination (-gaa) * eq0 + (-gab) * eq1 + ea * hD
    exact (mul_eq_zero.mp this).resolve_right hdt
  have I4 : eb - q0 * gab - q1 * gbb = 0 := by
    have : (eb - q0 * gab - q1 * gbb) * D = 0 := by
      linear_combination (-gab) * eq0 + (-gbb) * eq1 + eb * hD
    exact (mul_eq_zero.mp this).resolve_right hdt
  have hnew : new = (q0 + μ * p0) • Pa + (q1 + μ * p1) • Pb + (-μ) • Pd := by
    ext i; simp only [new, n, Pi.add_apply, Pi.sub_apply, Pi.smul_apply, smul_eq_mul]; ring
  have hba : bform M Pb Pa = gab := (hsym _ _).trans hab
  have had : bform M Pa Pd = gda := (hsym _ _).trans hda
  have hbd : bform M Pb Pd = gdb := (hsym _ _).trans hdb
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [hnew, bform_comb_left, haa, hba, hda]
    linear_combination (-1 : ℝ) * I3 + (-μ) * I1
  · rw [hnew, bform_comb_left, hab, hbb, hdb]
    linear_combination (-1 : ℝ) * I4 + (-μ) * I2
  · rw [hnew, bform_comb, haa, hab, had, hba, hbb, hbd, hda, hdb, hdd]
    linear_combination (-(μ ^ 2 * p0 + 2 * μ * q0)) * I1 + (-(μ ^ 2 * p1 + 2 * μ * q1)) * I2 +
      (-q0) * I3 + (-q1) * I4 + hμ
  · rw [hnew]
    simp only [rdet_eq, Pi.add_apply, Pi.smul_apply, smul_eq_mul]
    ring


/-! ### Crossings happen in strip order -/

lemma rdet_lin_right (A B X Y : Fin 3 → ℝ) (a b : ℝ) :
    rdet A B (a • X + b • Y) = a * rdet A B X + b * rdet A B Y := by
  simp only [rdet_eq, Pi.add_apply, Pi.smul_apply, smul_eq_mul]; ring

lemma rdet_self_left (A B : Fin 3 → ℝ) : rdet A B A = 0 := by simp only [rdet_eq]; ring
lemma rdet_self_right (A B : Fin 3 → ℝ) : rdet A B B = 0 := by simp only [rdet_eq]; ring
lemma rdet_swap12 (A B C : Fin 3 → ℝ) : rdet B A C = -rdet A B C := by simp only [rdet_eq]; ring
lemma rdet_cycle (A B C : Fin 3 → ℝ) : rdet A B C = rdet B C A := by simp only [rdet_eq]; ring

lemma rdet_ray (A B S E : Fin 3 → ℝ) (τ : ℝ) :
    rdet A B (S + τ • E) = rdet A B S + τ * rdet A B E := by
  simp only [rdet_eq, Pi.add_apply, Pi.smul_apply, smul_eq_mul]; ring

/-- One step of the ordering argument. Triangle `k` has the edge `A B` and
third vertex `D`; the next triangle has its new vertex `N` on the other side
of `A B` (`det(A,B,N) = -μ det(A,B,D)`) and the next edge is `N A` (so its
third vertex is `B`). If `S` lies on `D`'s side of `A B`, and the segment from
`S` toward `E` meets the open edge `A B` at parameter `τ` and the open edge
`N A` at `τ'`, then `τ < τ'` and `S` lies on `B`'s side of `N A`. -/
theorem order_step {A B D N S E : Fin 3 → ℝ} {μ τ τ' a1 b1 a2 b2 : ℝ}
    (hμ : 0 < μ) (hN : rdet A B N = -μ * rdet A B D)
    (hS : 0 < rdet A B S * rdet A B D)
    (hτ : 0 < τ) (ha1 : 0 < a1) (hb1 : 0 < b1) (h1 : S + τ • E = a1 • A + b1 • B)
    (ha2 : 0 < a2) (hb2 : 0 < b2) (h2 : S + τ' • E = a2 • N + b2 • A) :
    τ < τ' ∧ 0 < rdet N A S * rdet N A B := by
  set s := rdet A B S
  set dD := rdet A B D
  set e := rdet A B E
  have f1 : s + τ * e = 0 := by
    have := congrArg (rdet A B) h1
    rw [rdet_ray, rdet_lin_right, rdet_self_left, rdet_self_right] at this
    linarith
  have f2 : s + τ' * e = a2 * (-μ * dD) := by
    have := congrArg (rdet A B) h2
    rw [rdet_ray, rdet_lin_right, hN, rdet_self_left] at this
    linarith
  have hD0 : dD ≠ 0 := by intro h0; rw [h0, mul_zero] at hS; exact lt_irrefl _ hS
  have he : (τ' - τ) * e = -(a2 * μ * dD) := by linear_combination f2 - f1
  have hs : s = -τ * e := by linarith
  have key : (τ' - τ) * (s * dD) = τ * a2 * μ * dD ^ 2 := by
    rw [hs]; linear_combination (-τ * dD) * he
  have hpos : 0 < (τ' - τ) * (s * dD) := by
    rw [key]; have : 0 < dD ^ 2 := by positivity
    positivity
  have hlt : τ < τ' := by
    by_contra hle; push Not at hle; nlinarith
  refine ⟨hlt, ?_⟩
  set g1 := rdet N A S
  set gE := rdet N A E
  set gB := rdet N A B
  have g2 : g1 + τ' * gE = 0 := by
    have := congrArg (rdet N A) h2
    rw [rdet_ray, rdet_lin_right, rdet_self_left, rdet_self_right] at this
    linarith
  have gτ : g1 + τ * gE = b1 * gB := by
    have := congrArg (rdet N A) h1
    rw [rdet_ray, rdet_lin_right, rdet_self_right] at this
    linarith
  have hgB : gB ≠ 0 := by
    have : gB = rdet A B N := rdet_cycle N A B
    rw [this, hN]; exact mul_ne_zero (neg_ne_zero.mpr hμ.ne') hD0
  have hgE : gE ≠ 0 := by
    intro h0
    rw [h0, mul_zero, add_zero] at g2 gτ
    rw [g2] at gτ
    exact hgB ((mul_eq_zero.mp gτ.symm).resolve_left hb1.ne')
  have hτ' : 0 < τ' := lt_trans hτ hlt
  have hprod : g1 * (b1 * gB) = τ' * (τ' - τ) * gE ^ 2 := by
    have hg1 : g1 = -τ' * gE := by linarith
    rw [← gτ, hg1]; ring
  have : 0 < g1 * (b1 * gB) := by
    rw [hprod]; have : 0 < gE ^ 2 := by positivity
    have : 0 < τ' - τ := by linarith
    positivity
  have : 0 < b1 * (g1 * gB) := by nlinarith
  exact pos_of_mul_pos_right this hb1.le

/-- The same with the next edge `N B` (third vertex `A`). -/
theorem order_step' {A B D N S E : Fin 3 → ℝ} {μ τ τ' a1 b1 a2 b2 : ℝ}
    (hμ : 0 < μ) (hN : rdet A B N = -μ * rdet A B D)
    (hS : 0 < rdet A B S * rdet A B D)
    (hτ : 0 < τ) (ha1 : 0 < a1) (hb1 : 0 < b1) (h1 : S + τ • E = a1 • A + b1 • B)
    (ha2 : 0 < a2) (hb2 : 0 < b2) (h2 : S + τ' • E = a2 • N + b2 • B) :
    τ < τ' ∧ 0 < rdet N B S * rdet N B A := by
  refine order_step (A := B) (B := A) (D := D) hμ ?_ ?_ hτ hb1 ha1 ?_ ha2 hb2 h2
  · rw [rdet_swap12, hN, rdet_swap12]; ring
  · rw [rdet_swap12, rdet_swap12 A B D]; nlinarith
  · rw [h1, add_comm]


/-- `order_step` for any position of the next edge: it joins the new vertex
`N` to one end of `A B` (in either order), and its third vertex `Z` is the
other end. -/
theorem order_next {A B D N S E X Y Z : Fin 3 → ℝ} {μ τ τ' a1 b1 a2 b2 : ℝ}
    (hμ : 0 < μ) (hN : rdet A B N = -μ * rdet A B D)
    (hS : 0 < rdet A B S * rdet A B D)
    (hτ : 0 < τ) (ha1 : 0 < a1) (hb1 : 0 < b1) (h1 : S + τ • E = a1 • A + b1 • B)
    (hcase : (X = A ∧ Y = N ∧ Z = B) ∨ (X = N ∧ Y = A ∧ Z = B) ∨
      (X = B ∧ Y = N ∧ Z = A) ∨ (X = N ∧ Y = B ∧ Z = A))
    (ha2 : 0 < a2) (hb2 : 0 < b2) (h2 : S + τ' • E = a2 • X + b2 • Y) :
    τ < τ' ∧ 0 < rdet X Y S * rdet X Y Z := by
  rcases hcase with ⟨hX, hY, hZ⟩ | ⟨hX, hY, hZ⟩ | ⟨hX, hY, hZ⟩ | ⟨hX, hY, hZ⟩ <;>
    rw [hX, hY] at h2 <;> rw [hX, hY, hZ]
  · obtain ⟨hlt, hpos⟩ := order_step hμ hN hS hτ ha1 hb1 h1 hb2 ha2 (h2.trans (add_comm _ _))
    refine ⟨hlt, ?_⟩
    rw [rdet_swap12 N A S, rdet_swap12 N A B, neg_mul_neg]; exact hpos
  · exact order_step hμ hN hS hτ ha1 hb1 h1 ha2 hb2 h2
  · obtain ⟨hlt, hpos⟩ := order_step' hμ hN hS hτ ha1 hb1 h1 hb2 ha2 (h2.trans (add_comm _ _))
    refine ⟨hlt, ?_⟩
    rw [rdet_swap12 N B S, rdet_swap12 N B A, neg_mul_neg]; exact hpos
  · exact order_step' hμ hN hS hτ ha1 hb1 h1 ha2 hb2 h2

end VLink
