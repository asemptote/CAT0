import VLink.MQ

/-!
# Nested square roots (layer L1, mirrors `Tower` in `vlink/field.py`)

For generic edge lengths the sine of an angle is `√r` with `r` a multiquadratic
number that is not rational. `Tw` adjoins such square roots one at a time:
`ext r A B` is `A + B √r`. As for `MQ`, no normal form is assumed; the order in
which radicands are nested (`outer`) only affects efficiency.

Multiplication uses `√r * √r = r`, which needs `r ≥ 0`: that is the
well-formedness condition `WF`, established when a square root is created
(`sqrtField`) and preserved by the operations.
-/

namespace VLink

/-! ### Square roots of rationals -/

/-- `n = k² s`, by trial division (any answer is fine: `sqrtRat` checks it). -/
def sqfree (n : ℕ) : ℕ × ℕ := go n n 1 1 2
where
  go : ℕ → ℕ → ℕ → ℕ → ℕ → ℕ × ℕ
    | 0, n, k, s, _ => (k, s * n)
    | fuel + 1, n, k, s, p =>
      if n ≤ 1 then (k, s)
      else if p * p > n then (k, s * n)
      else if n % (p * p) = 0 then go fuel (n / (p * p)) (k * p) s p
      else if n % p = 0 then go fuel (n / p) k (s * p) (p + 1)
      else go fuel n k s (p + 1)

/-- `√q` for a rational `q ≥ 0`: `√(n/d) = √(n d) / d = (k/d) √s`. -/
def MQ.sqrtRat (q : ℚ) : MQ :=
  let n := q.num.toNat
  let d := q.den
  let ks := sqfree (n * d)
  if ks.1 * ks.1 * ks.2 = n * d then [(ks.2, (ks.1 : ℚ) / d)] else [(n * d, 1 / (d : ℚ))]

lemma sqrt_rat_eq (q : ℚ) (hq : 0 ≤ q) :
    Real.sqrt (q : ℝ) = Real.sqrt ((q.num.toNat * q.den : ℕ) : ℝ) / q.den := by
  have hd : (0 : ℝ) < q.den := by exact_mod_cast q.den_pos
  have hn : ((q.num.toNat : ℕ) : ℝ) = (q.num : ℝ) := by
    exact_mod_cast Int.toNat_of_nonneg (Rat.num_nonneg.mpr hq)
  have hqd : (q : ℝ) = (q.num.toNat * q.den : ℕ) / ((q.den : ℝ) * q.den) := by
    push_cast
    rw [hn, Rat.cast_def]
    field_simp
  rw [hqd, Real.sqrt_div' _ (by positivity), Real.sqrt_mul_self hd.le]

lemma MQ.val_sqrtRat (q : ℚ) (hq : 0 ≤ q) : MQ.val (MQ.sqrtRat q) = Real.sqrt q := by
  rw [sqrt_rat_eq q hq]
  unfold MQ.sqrtRat
  dsimp only
  generalize sqfree (q.num.toNat * q.den) = ks
  split_ifs with h
  · have e : ((q.num.toNat * q.den : ℕ) : ℝ) = (ks.1 : ℝ) ^ 2 * ks.2 := by
      rw [← h]; push_cast; ring
    rw [MQ.val_cons, MQ.val_nil, e, Real.sqrt_mul (by positivity),
      Real.sqrt_sq (Nat.cast_nonneg _)]
    push_cast; ring
  · rw [MQ.val_cons, MQ.val_nil]; push_cast; ring

/-! ### Rational enclosures of square roots -/

/-- Lower bound of `√q` (`0` if `q ≤ 0`). -/
def qsqrtLo (q : ℚ) : ℚ := if q ≤ 0 then 0 else MQ.sqrtLo (q.num.toNat * q.den) / q.den

/-- Upper bound of `√q` (`0` if `q ≤ 0`). -/
def qsqrtHi (q : ℚ) : ℚ := if q ≤ 0 then 0 else MQ.sqrtHi (q.num.toNat * q.den) / q.den

lemma qsqrtLo_nonneg (q : ℚ) : 0 ≤ qsqrtLo q := by
  unfold qsqrtLo MQ.sqrtLo
  split_ifs <;> positivity

lemma qsqrtLo_le (q : ℚ) : ((qsqrtLo q : ℚ) : ℝ) ≤ Real.sqrt q := by
  unfold qsqrtLo
  split_ifs with h
  · simp [Real.sqrt_nonneg]
  · rw [sqrt_rat_eq q (le_of_lt (not_le.mp h))]
    have := MQ.sqrtLo_le (q.num.toNat * q.den)
    push_cast at this ⊢
    exact div_le_div_of_nonneg_right this (Nat.cast_nonneg _)

lemma le_qsqrtHi (q : ℚ) : Real.sqrt q ≤ ((qsqrtHi q : ℚ) : ℝ) := by
  unfold qsqrtHi
  split_ifs with h
  · rw [Real.sqrt_eq_zero'.mpr (by exact_mod_cast h)]; simp
  · rw [sqrt_rat_eq q (le_of_lt (not_le.mp h))]
    have := MQ.le_sqrtHi (q.num.toNat * q.den)
    push_cast at this ⊢
    exact div_le_div_of_nonneg_right this (Nat.cast_nonneg _)

/-- Interval product, second factor non-negative. -/
lemma mul_mem_iv {x y l1 h1 l2 h2 : ℝ} (hx : l1 ≤ x ∧ x ≤ h1) (hy : l2 ≤ y ∧ y ≤ h2)
    (hl2 : 0 ≤ l2) : min (l1 * l2) (l1 * h2) ≤ x * y ∧ x * y ≤ max (h1 * l2) (h1 * h2) := by
  have y0 : 0 ≤ y := hl2.trans hy.1
  constructor
  · rcases le_or_gt 0 l1 with h | h
    · exact (min_le_left _ _).trans (by nlinarith)
    · exact (min_le_right _ _).trans (by nlinarith)
  · rcases le_or_gt 0 h1 with h | h
    · exact (le_max_right _ _).trans' (by nlinarith)
    · exact (le_max_left _ _).trans' (by nlinarith)

/-! ### Towers -/

/-- `base x` is `x`; `ext r A B` is `A + B √r`. -/
inductive Tw where
  | base : MQ → Tw
  | ext : MQ → Tw → Tw → Tw
  deriving Inhabited

namespace Tw

noncomputable def val : Tw → ℝ
  | base x => x.val
  | ext r A B => val A + val B * Real.sqrt r.val

/-- Every radicand is non-negative. -/
def WF : Tw → Prop
  | base _ => True
  | ext r A B => 0 ≤ r.val ∧ WF A ∧ WF B

/-- Lexicographic order on radicands; decides which square root is outermost. -/
def cmpR : MQ → MQ → Ordering
  | [], [] => .eq
  | [], _ :: _ => .lt
  | _ :: _, [] => .gt
  | p :: l, q :: m =>
    match compare p.1 q.1 with
    | .eq => match compare p.2 q.2 with
      | .eq => cmpR l m
      | o => o
    | o => o

def outer (r s : MQ) : Bool := cmpR r s == .gt

def add : Tw → Tw → Tw
  | base x, base y => base (MQ.add x y)
  | base x, ext s C D => ext s (add (base x) C) D
  | ext r A B, base y => ext r (add A (base y)) B
  | ext r A B, ext s C D =>
    if r = s then ext r (add A C) (add B D)
    else if outer r s then ext r (add A (ext s C D)) B
    else ext s (add (ext r A B) C) D

lemma val_add (x y : Tw) : val (add x y) = val x + val y := by
  fun_induction add x y <;> simp_all [add, val, MQ.val_add] <;> ring

lemma WF_add {x y : Tw} (hx : WF x) (hy : WF y) : WF (add x y) := by
  fun_induction add x y <;> simp_all [add, WF]

/-- Multiplication by a multiquadratic number. -/
def mulBase (m : MQ) : Tw → Tw
  | base x => base (MQ.mul m x)
  | ext r A B => ext r (mulBase m A) (mulBase m B)

lemma val_mulBase (m : MQ) (x : Tw) : val (mulBase m x) = m.val * val x := by
  induction x with
  | base x => simp [mulBase, val, MQ.val_mul]
  | ext r A B iA iB => simp only [mulBase, val, iA, iB]; ring

lemma WF_mulBase (m : MQ) {x : Tw} (hx : WF x) : WF (mulBase m x) := by
  induction x with
  | base x => trivial
  | ext r A B iA iB => exact ⟨hx.1, iA hx.2.1, iB hx.2.2⟩

def mul : Tw → Tw → Tw
  | base x, y => mulBase x y
  | ext r A B, base y => ext r (mulBase y A) (mulBase y B)
  | ext r A B, ext s C D =>
    if r = s then ext r (add (mul A C) (mulBase r (mul B D))) (add (mul A D) (mul B C))
    else if outer r s then ext r (mul A (ext s C D)) (mul B (ext s C D))
    else ext s (mul (ext r A B) C) (mul (ext r A B) D)

lemma WF_mul {x y : Tw} (hx : WF x) (hy : WF y) : WF (mul x y) := by
  fun_induction mul x y <;>
    simp_all [mul, WF, WF_add, WF_mulBase]

theorem val_mul : ∀ {x y : Tw}, WF x → WF y → val (mul x y) = val x * val y
  | base x, y, _, _ => by simp [mul, val, val_mulBase]
  | ext r A B, base y, _, _ => by simp only [mul, val, val_mulBase]; ring
  | ext r A B, ext s C D, hx, hy => by
    have hx' := hx
    have hy' := hy
    obtain ⟨hr, hA, hB⟩ := hx
    obtain ⟨hs, hC, hD⟩ := hy
    by_cases h : r = s
    · subst h
      have iAC := val_mul hA hC
      have iBD := val_mul hB hD
      have iAD := val_mul hA hD
      have iBC := val_mul hB hC
      simp only [mul, ↓reduceIte, ite_true, val, val_add, val_mulBase, iAC, iBD, iAD, iBC]
      have := Real.mul_self_sqrt hr
      linear_combination (-(val B * val D)) * this
    · by_cases ho : outer r s = true
      · have iA := val_mul hA hy'
        have iB := val_mul hB hy'
        simp only [mul, if_neg h, ho, if_true, val, iA, iB]
        ring
      · have iC := val_mul hx' hC
        have iD := val_mul hx' hD
        simp only [mul, if_neg h, ho, Bool.false_eq_true, if_false, val, iC, iD]
        ring

def sub (x y : Tw) : Tw := add x (mulBase (MQ.ofRat (-1)) y)

lemma val_sub (x y : Tw) : val (sub x y) = val x - val y := by
  simp [sub, val_add, val_mulBase]; ring

lemma WF_sub {x y : Tw} (hx : WF x) (hy : WF y) : WF (sub x y) :=
  WF_add hx (WF_mulBase _ hy)

/-! ### Sign -/

/-- A rational interval containing `val x`. -/
def iv : Tw → ℚ × ℚ
  | base x => MQ.iv x
  | ext r A B =>
    let a := iv A
    let b := iv B
    let lo := qsqrtLo (MQ.iv r).1
    let hi := qsqrtHi (MQ.iv r).2
    (a.1 + min (b.1 * lo) (b.1 * hi), a.2 + max (b.2 * lo) (b.2 * hi))

lemma iv_sound (x : Tw) : (((iv x).1 : ℚ) : ℝ) ≤ val x ∧ val x ≤ (((iv x).2 : ℚ) : ℝ) := by
  induction x with
  | base x => exact MQ.iv_sound x
  | ext r A B iA iB =>
    have hr := MQ.iv_sound r
    have hs : ((qsqrtLo (MQ.iv r).1 : ℚ) : ℝ) ≤ Real.sqrt r.val ∧
        Real.sqrt r.val ≤ ((qsqrtHi (MQ.iv r).2 : ℚ) : ℝ) :=
      ⟨(qsqrtLo_le _).trans (Real.sqrt_le_sqrt hr.1),
       (Real.sqrt_le_sqrt hr.2).trans (le_qsqrtHi _)⟩
    have hm := mul_mem_iv iB hs (by exact_mod_cast qsqrtLo_nonneg _)
    simp only [iv, val]
    push_cast
    constructor <;> linarith [iA.1, iA.2, hm.1, hm.2]

/-- Fuel for the multiquadratic signs called from here. -/
def mqFuel : ℕ := 64

/-- Exact sign of `val x` (for `WF x`), or `none` when out of fuel. -/
def sign : ℕ → Tw → Option ℤ
  | 0, _ => none
  | _ + 1, base x => MQ.sign mqFuel x
  | fuel + 1, ext r A B =>
    if 0 < (iv (ext r A B)).1 then some 1
    else if (iv (ext r A B)).2 < 0 then some (-1)
    else match MQ.sign mqFuel r with
      | none => none
      | some sr =>
        if sr ≤ 0 then sign fuel A      -- r ≤ 0: √r = 0
        else match sign fuel A, sign fuel B with
          | some sa, some sb =>
            if sb = 0 then some sa
            else if sa = 0 then some sb
            else if sa = sb then some sa
            else (sign fuel (sub (mul A A) (mulBase r (mul B B)))).map (sa * ·)
          | _, _ => none

theorem sign_sound : ∀ (fuel : ℕ) (x : Tw) (s : ℤ), WF x → sign fuel x = some s →
    s = rsign (val x)
  | 0, x, s, _, h => by simp [sign] at h
  | _ + 1, base x, s, _, h => MQ.sign_sound _ x s h
  | fuel + 1, ext r A B, s, hw, h => by
    have ih := sign_sound fuel
    have hiv := iv_sound (ext r A B)
    obtain ⟨hr0, hA, hB⟩ := hw
    unfold sign at h
    by_cases h1 : 0 < (iv (ext r A B)).1
    · rw [if_pos h1] at h
      cases h
      exact (rsign_of_pos (lt_of_lt_of_le (by exact_mod_cast h1) hiv.1)).symm
    rw [if_neg h1] at h
    by_cases h2 : (iv (ext r A B)).2 < 0
    · rw [if_pos h2] at h
      cases h
      exact (rsign_of_neg (lt_of_le_of_lt hiv.2 (by exact_mod_cast h2))).symm
    rw [if_neg h2] at h
    cases hsr : MQ.sign mqFuel r with
    | none => simp [hsr] at h
    | some sr =>
    simp only [hsr] at h
    have er := MQ.sign_sound _ r sr hsr
    simp only [val]
    by_cases c0 : sr ≤ 0
    · rw [if_pos c0] at h
      have : r.val ≤ 0 := by
        rcases rsign_cases r.val with ⟨_, e⟩ | ⟨h', _⟩ | ⟨h', _⟩ <;> [omega; exact h'.le; exact h'.le]
      rw [Real.sqrt_eq_zero'.mpr this, mul_zero, add_zero]
      exact ih A s hA h
    rw [if_neg c0] at h
    have hrpos : 0 < r.val := by
      rcases rsign_cases r.val with ⟨h', _⟩ | ⟨_, e⟩ | ⟨_, e⟩ <;> [exact h'; omega; omega]
    have hsq := Real.sqrt_nonneg r.val
    cases ha : sign fuel A with
    | none => simp [ha] at h
    | some sa =>
    cases hb : sign fuel B with
    | none => simp [ha, hb] at h
    | some sb =>
    simp only [ha, hb] at h
    have ea := ih A sa hA ha
    have eb := ih B sb hB hb
    by_cases c1 : sb = 0
    · rw [if_pos c1] at h
      cases h
      rw [rsign_eq_zero.mp (eb.symm.trans c1), zero_mul, add_zero, ea]
    rw [if_neg c1] at h
    by_cases c2 : sa = 0
    · rw [if_pos c2] at h
      cases h
      rw [rsign_eq_zero.mp (ea.symm.trans c2), zero_add, eb,
        rsign_mul_pos (Real.sqrt_pos.mpr hrpos)]
    rw [if_neg c2] at h
    by_cases c3 : sa = sb
    · rw [if_pos c3] at h
      cases h
      rw [ea]
      rcases rsign_cases (val A) with ⟨hA', e⟩ | ⟨hA', e⟩ | ⟨hA', e⟩
      · have hB' : 0 < val B := by
          rcases rsign_cases (val B) with ⟨hB', _⟩ | ⟨_, f⟩ | ⟨_, f⟩
          · exact hB'
          all_goals omega
        rw [e, rsign_of_pos (by nlinarith)]
      · have hB' : val B < 0 := by
          rcases rsign_cases (val B) with ⟨_, f⟩ | ⟨hB', _⟩ | ⟨_, f⟩
          · omega
          · exact hB'
          · omega
        rw [e, rsign_of_neg (by nlinarith)]
      · omega
    rw [if_neg c3] at h
    obtain ⟨sd, hd, rfl⟩ := Option.map_eq_some_iff.mp h
    have ed := ih _ sd (WF_sub (WF_mul hA hA) (WF_mulBase _ (WF_mul hB hB))) hd
    rw [val_sub, val_mulBase, val_mul hA hA, val_mul hB hB] at ed
    rw [ea, ed, rsign_split hrpos
      (fun h0 => c2 (by rw [ea, h0, rsign_zero]))
      (fun h0 => c1 (by rw [eb, h0, rsign_zero]))
      (by rw [← ea, ← eb]; exact c3)]
    congr 2
    ring

/-- `√b2` for `b2 ≥ 0`: in `MQ` when `b2` is rational, a new square root otherwise. -/
def sqrtField (b2 : MQ) : Option Tw :=
  let b := MQ.norm b2
  if b.all (fun p => p.1 ≤ 1) then
    let q := MQ.ratPart b
    if q < 0 then none else some (base (MQ.sqrtRat q))
  else match MQ.sign mqFuel b with
    | some s => if 0 ≤ s then some (ext b (base []) (base (MQ.ofRat 1))) else none
    | none => none

lemma sqrtField_sound {b2 : MQ} {x : Tw} (h : sqrtField b2 = some x) :
    val x = Real.sqrt b2.val ∧ WF x := by
  unfold sqrtField at h
  simp only at h
  split_ifs at h with h1 h2
  · cases h
    have hall : ∀ p ∈ MQ.norm b2, p.1 ≤ 1 := by simpa using h1
    have e := MQ.val_ratPart _ hall
    rw [MQ.val_norm] at e
    refine ⟨?_, trivial⟩
    rw [val, MQ.val_sqrtRat _ (not_lt.mp h2), e]
  · split at h
    · rename_i s hs
      split_ifs at h with h3
      cases h
      have es := MQ.sign_sound _ _ s hs
      rw [MQ.val_norm] at es
      have hnn : 0 ≤ b2.val := by
        rcases rsign_cases b2.val with ⟨h', _⟩ | ⟨_, e⟩ | ⟨h', _⟩ <;> [exact h'.le; omega; exact h'.ge]
      refine ⟨?_, ⟨by rwa [MQ.val_norm], trivial, trivial⟩⟩
      simp [val, MQ.val_norm]
    · cases h

end Tw

end VLink
