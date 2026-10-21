import VLink.Tower

/-!
# Exact angle sums (layer L2, mirrors `AngleSum` in `vlink/kernel.py`)

An angle `θ ∈ [0, π]` is given by exact `a` and `b²` (`b ≥ 0`, `(a, b) ≠ 0`):
`θ = arg (a + i b)`. A sum of angles `Θ` is kept as `z = re + i im`, a positive
multiple of `exp (i Θ)`, together with the number `h` of half-turns:
`Θ ∈ [h π, (h + 1) π)`. Adding `θ` multiplies `z` by `a + i b` and moves `Θ`
into `[h π, (h + 2) π)`; which half is decided exactly from the signs of
`(-1)^h Im z` and `(-1)^h Re z`. Hence `Θ ≥ k π` iff `h ≥ k` (`atLeast_iff`).
-/

namespace VLink

open Real

/-- The angle `arg (a + i √b2) ∈ [0, π]`. -/
noncomputable def angleOf (a b2 : MQ) : ℝ := Complex.arg ⟨a.val, Real.sqrt b2.val⟩

lemma angleOf_spec (a b2 : MQ) (hne : a.val ≠ 0 ∨ 0 < b2.val) :
    ∃ σ > 0, a.val = σ * cos (angleOf a b2) ∧ Real.sqrt b2.val = σ * sin (angleOf a b2) ∧
      0 ≤ angleOf a b2 ∧ angleOf a b2 ≤ π := by
  set z : ℂ := ⟨a.val, Real.sqrt b2.val⟩ with hz
  have hz0 : z ≠ 0 := by
    intro h0
    have hre : a.val = 0 := by simpa [hz] using congrArg Complex.re h0
    have him : Real.sqrt b2.val = 0 := by simpa [hz] using congrArg Complex.im h0
    rcases hne with h | h
    · exact h hre
    · exact (Real.sqrt_pos.mpr h).ne' him
  refine ⟨‖z‖, norm_pos_iff.mpr hz0, ?_, ?_, ?_, Complex.arg_le_pi z⟩
  · have := Complex.norm_mul_cos_arg z
    simp only [angleOf, ← hz]
    rw [this]
  · have := Complex.norm_mul_sin_arg z
    simp only [angleOf, ← hz]
    rw [this]
  · exact Complex.arg_nonneg_iff.mpr (Real.sqrt_nonneg _)

/-- The angle as an arccosine: `θ = arccos (a / √(a² + b²))`. -/
lemma angleOf_eq_arccos (a b2 : MQ) (hb : 0 ≤ b2.val) (hne : a.val ≠ 0 ∨ 0 < b2.val) :
    angleOf a b2 = Real.arccos (a.val / Real.sqrt (a.val ^ 2 + b2.val)) := by
  obtain ⟨σ, hσ, ha, hb', h0, hπ⟩ := angleOf_spec a b2 hne
  have hs : Real.sqrt (a.val ^ 2 + b2.val) = σ := by
    have : a.val ^ 2 + b2.val = σ ^ 2 := by
      have e := Real.sq_sqrt hb
      rw [← e, ha, hb']
      have := sin_sq_add_cos_sq (angleOf a b2)
      linear_combination σ ^ 2 * this
    rw [this, Real.sqrt_sq hσ.le]
  rw [hs, ha, mul_div_cancel_left₀ _ hσ.ne', Real.arccos_cos h0 hπ]

structure AngleSum where
  re : Tw
  im : Tw
  h : ℕ

namespace AngleSum

/-- Fuel for the exact signs. -/
def fuel : ℕ := 64

def zero : AngleSum := ⟨.base (MQ.ofRat 1), .base [], 0⟩

/-- `st` represents the angle sum `Θ`. -/
def Rep (st : AngleSum) (Θ : ℝ) : Prop :=
  st.re.WF ∧ st.im.WF ∧ ∃ ρ > 0, st.re.val = ρ * cos Θ ∧ st.im.val = ρ * sin Θ ∧
    (st.h : ℝ) * π ≤ Θ ∧ Θ < ((st.h : ℝ) + 1) * π

lemma rep_zero : zero.Rep 0 :=
  ⟨trivial, trivial, 1, one_pos, by simp [zero, Tw.val], by simp [zero, Tw.val],
    by simp [zero], by simp [zero, pi_pos]⟩

/-- Add the angle `angleOf a b2`; `none` if the data is invalid or undecided. -/
def add (st : AngleSum) (a b2 : MQ) : Option AngleSum :=
  match Tw.sqrtField b2, MQ.sign fuel a, MQ.sign fuel b2 with
  | some b, some sa, some sb =>
    if sa = 0 ∧ sb ≤ 0 then none
    else
      let re := Tw.sub (Tw.mul st.re (.base a)) (Tw.mul st.im b)
      let im := Tw.add (Tw.mul st.re b) (Tw.mul st.im (.base a))
      match Tw.sign fuel im, Tw.sign fuel re with
      | some si, some sr =>
        let s : ℤ := if st.h % 2 = 0 then 1 else -1
        some ⟨re, im, if s * si < 0 ∨ (s * si = 0 ∧ s * sr < 0) then st.h + 1 else st.h⟩
      | _, _ => none
  | _, _, _ => none

lemma half_turn (ψ : ℝ) (h0 : 0 ≤ ψ) (h2 : ψ < 2 * π) :
    (sin ψ < 0 ∨ (sin ψ = 0 ∧ cos ψ < 0)) ↔ π ≤ ψ := by
  constructor
  · intro hc
    by_contra hlt
    push Not at hlt
    rcases h0.lt_or_eq with hpos | hz
    · have := sin_pos_of_pos_of_lt_pi hpos hlt
      rcases hc with hc | hc
      · linarith
      · linarith [hc.1]
    · subst hz
      norm_num at hc
  · intro hp
    rcases hp.lt_or_eq with hlt | heq
    · left
      have : sin (ψ - 2 * π) < 0 := sin_neg_of_neg_of_neg_pi_lt (by linarith) (by linarith)
      rwa [sin_sub_two_pi] at this
    · right
      rw [← heq]
      simp

lemma parity_sign (h : ℕ) : (((if h % 2 = 0 then 1 else -1 : ℤ)) : ℝ) = (-1) ^ h := by
  rcases Nat.mod_two_eq_zero_or_one h with e | e
  · rw [if_pos e, (Nat.even_iff.mpr e).neg_one_pow]; simp
  · rw [if_neg (by omega), (Nat.odd_iff.mpr e).neg_one_pow]; simp

lemma int_mul_rsign (s : ℤ) (hs : s = 1 ∨ s = -1) (u : ℝ) : s * rsign u = rsign ((s : ℝ) * u) := by
  rcases hs with rfl | rfl
  · simp
  · simp [rsign_neg]

theorem add_sound {st st' : AngleSum} {Θ : ℝ} {a b2 : MQ} (hr : st.Rep Θ)
    (h : st.add a b2 = some st') : st'.Rep (Θ + angleOf a b2) := by
  obtain ⟨wre, wim, ρ, hρ, hre, him, hlo, hhi⟩ := hr
  unfold add at h
  cases hb : Tw.sqrtField b2 with
  | none => simp [hb] at h
  | some b =>
  cases hsa : MQ.sign fuel a with
  | none => simp [hb, hsa] at h
  | some sa =>
  cases hsb : MQ.sign fuel b2 with
  | none => simp [hb, hsa, hsb] at h
  | some sb =>
  simp only [hb, hsa, hsb] at h
  obtain ⟨vb, wb⟩ := Tw.sqrtField_sound hb
  have ea := MQ.sign_sound _ _ _ hsa
  have eb := MQ.sign_sound _ _ _ hsb
  by_cases c : sa = 0 ∧ sb ≤ 0
  · rw [if_pos c] at h; cases h
  rw [if_neg c] at h
  have hne : a.val ≠ 0 ∨ 0 < b2.val := by
    by_contra hh
    push Not at hh
    apply c
    constructor
    · rw [ea, hh.1, rsign_zero]
    · rw [eb]; rcases rsign_cases b2.val with ⟨h', _⟩ | ⟨_, e⟩ | ⟨_, e⟩
      · linarith [hh.2]
      · omega
      · omega
  obtain ⟨σ, hσ, ha, hbσ, h0, hπ⟩ := angleOf_spec a b2 hne
  set θ := angleOf a b2
  set re := Tw.sub (Tw.mul st.re (.base a)) (Tw.mul st.im b)
  set im := Tw.add (Tw.mul st.re b) (Tw.mul st.im (.base a))
  have wre' : re.WF := Tw.WF_sub (Tw.WF_mul wre trivial) (Tw.WF_mul wim wb)
  have wim' : im.WF := Tw.WF_add (Tw.WF_mul wre wb) (Tw.WF_mul wim trivial)
  have vre : re.val = ρ * σ * cos (Θ + θ) := by
    simp only [re, Tw.val_sub, Tw.val_mul wre (show Tw.WF (.base a) from trivial), Tw.val_mul wim wb, vb, Tw.val, hre, him,
      ha, hbσ, cos_add]
    ring
  have vim : im.val = ρ * σ * sin (Θ + θ) := by
    simp only [im, Tw.val_add, Tw.val_mul wre wb, Tw.val_mul wim (show Tw.WF (.base a) from trivial), vb, Tw.val, hre, him,
      ha, hbσ, sin_add]
    ring
  cases hsi : Tw.sign fuel im with
  | none => simp [hsi] at h
  | some si =>
  cases hsr : Tw.sign fuel re with
  | none => simp [hsi, hsr] at h
  | some sr =>
  simp only [hsi, hsr, Option.some.injEq] at h
  subst h
  have esi := Tw.sign_sound _ _ _ wim' hsi
  have esr := Tw.sign_sound _ _ _ wre' hsr
  refine ⟨wre', wim', ρ * σ, mul_pos hρ hσ, vre, vim, ?_⟩
  -- the half-turn update
  set s : ℤ := if st.h % 2 = 0 then 1 else -1 with hsdef
  have hs : s = 1 ∨ s = -1 := by rw [hsdef]; split_ifs <;> simp
  set ψ := Θ + θ - st.h * π
  have hψ0 : 0 ≤ ψ := by simp only [ψ]; linarith
  have hψ2 : ψ < 2 * π := by simp only [ψ]; linarith
  have hsin : (s : ℝ) * sin (Θ + θ) = sin ψ := by
    have e : Θ + θ = ψ + st.h * π := by simp only [ψ]; ring
    rw [e, sin_add_nat_mul_pi, hsdef, parity_sign, ← mul_assoc, ← pow_add, ← two_mul, pow_mul]
    simp
  have hcos : (s : ℝ) * cos (Θ + θ) = cos ψ := by
    have e : Θ + θ = ψ + st.h * π := by simp only [ψ]; ring
    rw [e, cos_add_nat_mul_pi, hsdef, parity_sign, ← mul_assoc, ← pow_add, ← two_mul, pow_mul]
    simp
  have k : 0 < ρ * σ := mul_pos hρ hσ
  have csi : s * si = rsign (sin ψ) := by
    rw [esi, vim, int_mul_rsign s hs, ← hsin, mul_comm (ρ * σ), ← mul_assoc,
      rsign_mul_pos k]
  have csr : s * sr = rsign (cos ψ) := by
    rw [esr, vre, int_mul_rsign s hs, ← hcos, mul_comm (ρ * σ), ← mul_assoc,
      rsign_mul_pos k]
  have cond : (s * si < 0 ∨ (s * si = 0 ∧ s * sr < 0)) ↔ π ≤ ψ := by
    rw [csi, csr, rsign_lt_zero, rsign_eq_zero, rsign_lt_zero]
    exact half_turn ψ hψ0 hψ2
  split_ifs with hc
  · have := cond.mp hc
    push_cast
    constructor <;> [skip; skip] <;> simp only [ψ] at this hψ2 <;> linarith
  · have := mt cond.mpr hc
    push Not at this
    constructor <;> simp only [ψ] at this hψ0 <;> linarith

/-- Add a list of angles `(a, b²)`. -/
def addList (st : AngleSum) : List (MQ × MQ) → Option AngleSum
  | [] => some st
  | p :: l => match st.add p.1 p.2 with
    | some st' => addList st' l
    | none => none

theorem addList_sound {st st' : AngleSum} {Θ : ℝ} {l : List (MQ × MQ)} (hr : st.Rep Θ)
    (h : st.addList l = some st') :
    st'.Rep (Θ + (l.map fun p => angleOf p.1 p.2).sum) := by
  induction l generalizing st Θ with
  | nil => simp only [addList, Option.some.injEq] at h; subst h; simpa using hr
  | cons p l ih =>
    simp only [addList] at h
    cases ha : st.add p.1 p.2 with
    | none => simp [ha] at h
    | some st1 =>
      simp only [ha] at h
      have := ih (add_sound hr ha) h
      simpa [add_assoc] using this

theorem atLeast_iff {st : AngleSum} {Θ : ℝ} (hr : st.Rep Θ) (k : ℕ) :
    k ≤ st.h ↔ (k : ℝ) * π ≤ Θ := by
  obtain ⟨-, -, -, -, -, -, hlo, hhi⟩ := hr
  constructor
  · intro hk
    have : (k : ℝ) ≤ st.h := by exact_mod_cast hk
    nlinarith [pi_pos]
  · intro hk
    by_contra hlt
    push Not at hlt
    have : (st.h : ℝ) + 1 ≤ k := by exact_mod_cast hlt
    nlinarith [pi_pos]

/-- `Θ₁ ≤ Θ₂`, decided exactly (`angle_leq`). -/
def leq (x y : AngleSum) : Option Bool :=
  if x.h ≠ y.h then some (decide (x.h < y.h))
  else (Tw.sign fuel (Tw.sub (Tw.mul x.re y.im) (Tw.mul x.im y.re))).map (fun s => decide (0 ≤ s))

theorem leq_sound {x y : AngleSum} {Θ₁ Θ₂ : ℝ} (hx : x.Rep Θ₁) (hy : y.Rep Θ₂) {b : Bool}
    (h : x.leq y = some b) : (b = true ↔ Θ₁ ≤ Θ₂) := by
  obtain ⟨wx1, wx2, ρ, hρ, hre1, him1, hlo1, hhi1⟩ := hx
  obtain ⟨wy1, wy2, κ, hκ, hre2, him2, hlo2, hhi2⟩ := hy
  unfold leq at h
  split_ifs at h with hne
  · simp only [Option.some.injEq] at h
    subst h
    rcases Nat.lt_or_gt_of_ne hne with hl | hl
    · have : (x.h : ℝ) + 1 ≤ y.h := by exact_mod_cast hl
      simp only [hl, decide_true, true_iff]
      nlinarith [pi_pos]
    · have : (y.h : ℝ) + 1 ≤ x.h := by exact_mod_cast hl
      simp only [Nat.not_lt.mpr hl.le, decide_false, Bool.false_eq_true, false_iff, not_le]
      nlinarith [pi_pos]
  · push Not at hne
    obtain ⟨s, hs, rfl⟩ := Option.map_eq_some_iff.mp h
    have es := Tw.sign_sound _ _ _ (Tw.WF_sub (Tw.WF_mul wx1 wy2) (Tw.WF_mul wx2 wy1)) hs
    rw [Tw.val_sub, Tw.val_mul wx1 wy2, Tw.val_mul wx2 wy1, hre1, him1, hre2, him2] at es
    have e : ρ * cos Θ₁ * (κ * sin Θ₂) - ρ * sin Θ₁ * (κ * cos Θ₂) =
        sin (Θ₂ - Θ₁) * (ρ * κ) := by rw [sin_sub]; ring
    rw [e, rsign_mul_pos (mul_pos hρ hκ)] at es
    have hne' : (x.h : ℝ) = y.h := by exact_mod_cast hne
    simp only [decide_eq_true_eq, es, rsign_nonneg]
    constructor
    · intro hs0
      by_contra hlt
      push Not at hlt
      have := sin_neg_of_neg_of_neg_pi_lt (by linarith : Θ₂ - Θ₁ < 0) (by nlinarith [pi_pos])
      linarith
    · intro hle
      exact sin_nonneg_of_nonneg_of_le_pi (by linarith) (by nlinarith [pi_pos])

end AngleSum

end VLink
