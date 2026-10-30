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

end VLink
