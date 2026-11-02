import Mathlib

/-!
# Multiquadratic numbers (layer L1, mirrors `vlink/field.py`)

An element `Σ c_s √s` is a list of pairs `(s, c_s)`. Nothing here assumes a
normal form, square-free keys, or the linear independence of square roots:
every operation is proved correct for the real value `val`, and `sign` is
proved sound (when it answers, the answer is the sign of `val`). Zero tests go
through `sign`.
-/

namespace VLink

/-- The sign of a real number, as an integer. -/
noncomputable def rsign (r : ℝ) : ℤ := if 0 < r then 1 else if r < 0 then -1 else 0

lemma rsign_cases (r : ℝ) :
    (0 < r ∧ rsign r = 1) ∨ (r < 0 ∧ rsign r = -1) ∨ (r = 0 ∧ rsign r = 0) := by
  unfold rsign
  rcases lt_trichotomy r 0 with h | h | h
  · right; left; exact ⟨h, by simp [h, not_lt.mpr h.le]⟩
  · right; right; exact ⟨h, by simp [h]⟩
  · left; exact ⟨h, by simp [h]⟩

lemma rsign_of_pos {r : ℝ} (h : 0 < r) : rsign r = 1 := by simp [rsign, h]

lemma rsign_of_neg {r : ℝ} (h : r < 0) : rsign r = -1 := by
  simp [rsign, h, not_lt.mpr h.le]

@[simp] lemma rsign_zero : rsign 0 = 0 := by simp [rsign]

lemma rsign_eq_zero {r : ℝ} : rsign r = 0 ↔ r = 0 := by
  rcases rsign_cases r with ⟨h, e⟩ | ⟨h, e⟩ | ⟨h, e⟩ <;> rw [e]
  · simp [h.ne']
  · simp [h.ne]
  · simp [h]

lemma rsign_mul_pos {x y : ℝ} (hy : 0 < y) : rsign (x * y) = rsign x := by
  rcases rsign_cases x with ⟨h, e⟩ | ⟨h, e⟩ | ⟨h, e⟩ <;> rw [e]
  · exact rsign_of_pos (mul_pos h hy)
  · exact rsign_of_neg (mul_neg_of_neg_of_pos h hy)
  · simp [h]

lemma rsign_neg (r : ℝ) : rsign (-r) = -rsign r := by
  rcases rsign_cases r with ⟨h, e⟩ | ⟨h, e⟩ | ⟨h, e⟩ <;> rw [e]
  · exact rsign_of_neg (neg_lt_zero.mpr h)
  · rw [rsign_of_pos (neg_pos.mpr h)]; rfl
  · simp [h]

lemma rsign_lt_zero {r : ℝ} : rsign r < 0 ↔ r < 0 := by
  rcases rsign_cases r with ⟨h, e⟩ | ⟨h, e⟩ | ⟨h, e⟩ <;> rw [e]
  · simp [not_lt.mpr h.le]
  · simp [h]
  · simp [h]

lemma rsign_nonneg {r : ℝ} : 0 ≤ rsign r ↔ 0 ≤ r := by
  rcases rsign_cases r with ⟨h, e⟩ | ⟨h, e⟩ | ⟨h, e⟩ <;> rw [e]
  · simp [h.le]
  · simp [h]
  · simp [h]

/-- The heart of the sign algorithm: if `A` and `B` are nonzero with opposite
signs, `A + B √P` has the sign of `A` times the sign of `A² - B² P`. -/
lemma rsign_split {A B P : ℝ} (hP : 0 < P) (hA : A ≠ 0) (hB : B ≠ 0)
    (hopp : rsign A ≠ rsign B) :
    rsign (A + B * Real.sqrt P) = rsign A * rsign (A ^ 2 - B ^ 2 * P) := by
  have hq : 0 < Real.sqrt P := Real.sqrt_pos.mpr hP
  have hqq : Real.sqrt P * Real.sqrt P = P := Real.mul_self_sqrt hP.le
  have key : (A + B * Real.sqrt P) * (A - B * Real.sqrt P) = A ^ 2 - B ^ 2 * P := by
    linear_combination (-B ^ 2) * hqq
  rcases rsign_cases A with ⟨ha, ea⟩ | ⟨ha, ea⟩ | ⟨ha, ea⟩
  · -- A > 0, so B < 0 and A - B √P > 0
    rcases rsign_cases B with ⟨hb, eb⟩ | ⟨hb, eb⟩ | ⟨hb, eb⟩
    · exact absurd (ea.trans eb.symm) hopp
    · have hy : 0 < A - B * Real.sqrt P := by nlinarith
      rw [ea, ← key, rsign_mul_pos hy, one_mul]
    · exact absurd hb hB
  · -- A < 0, so B > 0 and A - B √P < 0
    rcases rsign_cases B with ⟨hb, eb⟩ | ⟨hb, eb⟩ | ⟨hb, eb⟩
    · have hy : 0 < -(A - B * Real.sqrt P) := by nlinarith
      have : A ^ 2 - B ^ 2 * P = -((A + B * Real.sqrt P) * -(A - B * Real.sqrt P)) := by
        rw [← key]; ring
      rw [ea, this, rsign_neg, rsign_mul_pos hy]; ring
    · exact absurd (ea.trans eb.symm) hopp
    · exact absurd hb hB
  · exact absurd ha hA

/-- The sign of a rational, computably. -/
def qsign (q : ℚ) : ℤ := if 0 < q then 1 else if q < 0 then -1 else 0

lemma qsign_eq (q : ℚ) : qsign q = rsign q := by
  unfold qsign rsign
  simp only [Rat.cast_pos, Rat.cast_lt_zero]

/-- `Σ c √s` as a list of `(s, c)`. -/
abbrev MQ := List (ℕ × ℚ)

namespace MQ

noncomputable def val : MQ → ℝ
  | [] => 0
  | p :: l => (p.2 : ℝ) * Real.sqrt p.1 + val l

@[simp] lemma val_nil : val [] = 0 := rfl

@[simp] lemma val_cons (p : ℕ × ℚ) (l : MQ) :
    val (p :: l) = (p.2 : ℝ) * Real.sqrt p.1 + val l := rfl

/-- The rational `q`. -/
def ofRat (q : ℚ) : MQ := if q = 0 then [] else [(1, q)]

@[simp] lemma val_ofRat (q : ℚ) : val (ofRat q) = q := by
  unfold ofRat; split_ifs with h <;> simp [h]

/-- Add the term `p` to `l`, merging it with a term of the same key. -/
def ins (p : ℕ × ℚ) : MQ → MQ
  | [] => if p.2 = 0 then [] else [p]
  | q :: l =>
    if q.1 = p.1 then (if q.2 + p.2 = 0 then l else (q.1, q.2 + p.2) :: l)
    else q :: ins p l

lemma val_ins (p : ℕ × ℚ) (l : MQ) :
    val (ins p l) = (p.2 : ℝ) * Real.sqrt p.1 + val l := by
  induction l with
  | nil => unfold ins; split_ifs with h <;> simp [h]
  | cons q l ih =>
    unfold ins
    split_ifs with h1 h2
    · have hc : (q.2 : ℝ) + p.2 = 0 := by exact_mod_cast h2
      simp only [val_cons]; rw [h1]
      linear_combination (-Real.sqrt (p.1 : ℝ)) * hc
    · simp only [val_cons]; rw [h1]; push_cast; ring
    · simp only [val_cons]; rw [ih]; ring

/-- Merge terms with equal keys and drop zero coefficients. -/
def norm (x : MQ) : MQ := x.foldr ins []

lemma val_norm (x : MQ) : val (norm x) = val x := by
  unfold norm
  induction x with
  | nil => rfl
  | cons p x ih => rw [List.foldr_cons, val_ins, ih, val_cons]

def add (x y : MQ) : MQ := y.foldr ins x

lemma val_add (x y : MQ) : val (add x y) = val x + val y := by
  unfold add
  induction y with
  | nil => simp
  | cons r y ih => rw [List.foldr_cons, val_ins, ih, val_cons]; ring

def smul (q : ℚ) (x : MQ) : MQ := if q = 0 then [] else x.map fun p => (p.1, q * p.2)

lemma val_smul (q : ℚ) (x : MQ) : val (smul q x) = q * val x := by
  unfold smul
  split_ifs with h
  · simp [h]
  · induction x with
    | nil => simp
    | cons p x ih => simp only [List.map_cons, val_cons, ih]; push_cast; ring

def neg (x : MQ) : MQ := x.map fun p => (p.1, -p.2)

lemma val_neg (x : MQ) : val (neg x) = -val x := by
  unfold neg
  induction x with
  | nil => simp
  | cons p x ih => simp only [List.map_cons, val_cons, ih]; push_cast; ring

def sub (x y : MQ) : MQ := add x (neg y)

lemma val_sub (x y : MQ) : val (sub x y) = val x - val y := by
  rw [sub, val_add, val_neg]; ring

/-- `(c₁ √s₁)(c₂ √s₂) = c₁ c₂ g √((s₁/g)(s₂/g))` with `g = gcd s₁ s₂`. -/
def mulTerm (p r : ℕ × ℚ) : ℕ × ℚ :=
  ((p.1 / Nat.gcd p.1 r.1) * (r.1 / Nat.gcd p.1 r.1), p.2 * r.2 * Nat.gcd p.1 r.1)

lemma val_mulTerm (p r : ℕ × ℚ) :
    ((mulTerm p r).2 : ℝ) * Real.sqrt (mulTerm p r).1 =
      ((p.2 : ℝ) * Real.sqrt p.1) * ((r.2 : ℝ) * Real.sqrt r.1) := by
  unfold mulTerm
  by_cases hg : Nat.gcd p.1 r.1 = 0
  · rw [Nat.gcd_eq_zero_iff] at hg
    simp [hg.1, hg.2]
  · obtain ⟨a, ha⟩ := Nat.gcd_dvd_left p.1 r.1
    obtain ⟨b, hb⟩ := Nat.gcd_dvd_right p.1 r.1
    have hpos : 0 < Nat.gcd p.1 r.1 := Nat.pos_of_ne_zero hg
    have hpa : p.1 / Nat.gcd p.1 r.1 = a :=
      Nat.div_eq_of_eq_mul_left hpos (by rw [mul_comm]; exact ha)
    have hrb : r.1 / Nat.gcd p.1 r.1 = b :=
      Nat.div_eq_of_eq_mul_left hpos (by rw [mul_comm]; exact hb)
    rw [hpa, hrb]
    generalize Nat.gcd p.1 r.1 = g at ha hb
    rw [ha, hb]
    push_cast
    have hg0 : (0 : ℝ) ≤ g := Nat.cast_nonneg g
    rw [Real.sqrt_mul (Nat.cast_nonneg a), Real.sqrt_mul hg0, Real.sqrt_mul hg0]
    have hgg : Real.sqrt (g : ℝ) * Real.sqrt g = g := Real.mul_self_sqrt hg0
    linear_combination (-((p.2 : ℝ) * r.2 * Real.sqrt a * Real.sqrt b)) * hgg

def mulRow (p : ℕ × ℚ) (y acc : MQ) : MQ := y.foldr (fun r a => ins (mulTerm p r) a) acc

lemma val_mulRow (p : ℕ × ℚ) (y acc : MQ) :
    val (mulRow p y acc) = ((p.2 : ℝ) * Real.sqrt p.1) * val y + val acc := by
  unfold mulRow
  induction y with
  | nil => simp
  | cons r y ih => rw [List.foldr_cons, val_ins, ih, val_mulTerm, val_cons]; ring

def mul (x y : MQ) : MQ := x.foldr (fun p acc => mulRow p y acc) []

lemma val_mul (x y : MQ) : val (mul x y) = val x * val y := by
  unfold mul
  induction x with
  | nil => simp
  | cons p x ih => rw [List.foldr_cons, val_mulRow, ih, val_cons]; ring

/-! ### Rigorous rational enclosures -/

/-- Working precision of the enclosures: `√s` is bracketed to within `1/D`. -/
def D : ℕ := 2 ^ 40

lemma D_pos : (0 : ℝ) < D := by unfold D; positivity

def sqrtLo (s : ℕ) : ℚ := (Nat.sqrt (s * D * D) : ℚ) / D

def sqrtHi (s : ℕ) : ℚ := ((Nat.sqrt (s * D * D) : ℚ) + 1) / D

lemma sqrtLo_le (s : ℕ) : ((sqrtLo s : ℚ) : ℝ) ≤ Real.sqrt s := by
  have hD := D_pos
  set m := Nat.sqrt (s * D * D)
  have hm : (m : ℝ) * m ≤ s * D * D := by exact_mod_cast Nat.sqrt_le (s * D * D)
  have h0 : (0 : ℝ) ≤ m / D := by positivity
  unfold sqrtLo
  push_cast
  calc (m : ℝ) / D = Real.sqrt (((m : ℝ) / D) ^ 2) := (Real.sqrt_sq h0).symm
    _ ≤ Real.sqrt s := by
      apply Real.sqrt_le_sqrt
      rw [div_pow, div_le_iff₀ (by positivity)]
      nlinarith

lemma le_sqrtHi (s : ℕ) : Real.sqrt s ≤ ((sqrtHi s : ℚ) : ℝ) := by
  have hD := D_pos
  set m := Nat.sqrt (s * D * D)
  have hm : (s : ℝ) * D * D < (m + 1) * (m + 1) := by
    exact_mod_cast Nat.lt_succ_sqrt (s * D * D)
  have h0 : (0 : ℝ) ≤ (m + 1) / D := by positivity
  unfold sqrtHi
  push_cast
  calc Real.sqrt s ≤ Real.sqrt ((((m : ℝ) + 1) / D) ^ 2) := by
        apply Real.sqrt_le_sqrt
        rw [div_pow, le_div_iff₀ (by positivity)]
        nlinarith
    _ = ((m : ℝ) + 1) / D := Real.sqrt_sq h0

def termIv (p : ℕ × ℚ) : ℚ × ℚ :=
  if 0 ≤ p.2 then (p.2 * sqrtLo p.1, p.2 * sqrtHi p.1)
  else (p.2 * sqrtHi p.1, p.2 * sqrtLo p.1)

lemma termIv_sound (p : ℕ × ℚ) :
    (((termIv p).1 : ℚ) : ℝ) ≤ (p.2 : ℝ) * Real.sqrt p.1 ∧
      (p.2 : ℝ) * Real.sqrt p.1 ≤ (((termIv p).2 : ℚ) : ℝ) := by
  have lo := sqrtLo_le p.1
  have hi := le_sqrtHi p.1
  unfold termIv
  split_ifs with h
  · have h' : (0 : ℝ) ≤ p.2 := by exact_mod_cast h
    push_cast
    exact ⟨mul_le_mul_of_nonneg_left lo h', mul_le_mul_of_nonneg_left hi h'⟩
  · have h' : (p.2 : ℝ) ≤ 0 := by exact_mod_cast (not_le.mp h).le
    push_cast
    exact ⟨mul_le_mul_of_nonpos_left hi h', mul_le_mul_of_nonpos_left lo h'⟩

/-- A rational interval containing `val x`. -/
def iv : MQ → ℚ × ℚ
  | [] => (0, 0)
  | p :: l => ((termIv p).1 + (iv l).1, (termIv p).2 + (iv l).2)

lemma iv_sound (x : MQ) : (((iv x).1 : ℚ) : ℝ) ≤ val x ∧ val x ≤ (((iv x).2 : ℚ) : ℝ) := by
  induction x with
  | nil => simp [iv]
  | cons p x ih =>
    have t := termIv_sound p
    simp only [iv, val_cons]
    push_cast
    constructor <;> linarith [t.1, t.2, ih.1, ih.2]

/-! ### Exact sign -/

/-- The rational part, for lists whose keys are all `0` or `1`. -/
def ratPart : MQ → ℚ
  | [] => 0
  | p :: l => (if p.1 = 1 then p.2 else 0) + ratPart l

lemma val_ratPart (x : MQ) (h : ∀ p ∈ x, p.1 ≤ 1) : val x = ratPart x := by
  induction x with
  | nil => simp [ratPart]
  | cons p x ih =>
    have hp := h p List.mem_cons_self
    rw [val_cons, ih (fun q hq => h q (List.mem_cons_of_mem p hq)), ratPart]
    rcases Nat.le_one_iff_eq_zero_or_eq_one.mp hp with e | e <;> simp [e]

/-- `x = A + B √k` with `A` the terms whose key `k` does not divide. -/
def split (k : ℕ) (x : MQ) : MQ × MQ :=
  (x.filter (fun p => decide (p.1 % k ≠ 0)),
   (x.filter (fun p => decide (p.1 % k = 0))).map (fun p => (p.1 / k, p.2)))

lemma val_split (k : ℕ) (x : MQ) :
    val x = val (split k x).1 + val (split k x).2 * Real.sqrt k := by
  unfold split
  dsimp only
  induction x with
  | nil => simp
  | cons p x ih =>
    by_cases h : p.1 % k = 0
    · have e : Real.sqrt (p.1 : ℝ) = Real.sqrt (k : ℝ) * Real.sqrt ((p.1 / k : ℕ) : ℝ) := by
        rw [← Real.sqrt_mul (Nat.cast_nonneg _), ← Nat.cast_mul,
          Nat.mul_div_cancel' (Nat.dvd_of_mod_eq_zero h)]
      simp only [List.filter_cons, h, ne_eq, not_true_eq_false, decide_false, decide_true,
        List.map_cons, val_cons, Bool.false_eq_true, ite_false, ite_true] at ih ⊢
      rw [ih, e]; ring
    · simp only [List.filter_cons, h, ne_eq, not_false_eq_true, decide_true, decide_false,
        Bool.false_eq_true, ite_false, ite_true, val_cons] at ih ⊢
      rw [ih]; ring

/-- Some `s > 1` that every key is divisible by or coprime to, by gcd
refinement (only used for efficiency: soundness holds for any `k ≥ 1`). -/
def refine (keys : List ℕ) : ℕ → ℕ → ℕ
  | 0, s => s
  | fuel + 1, s =>
    match keys.find? (fun k => 1 < Nat.gcd s k && Nat.gcd s k < s) with
    | some k => refine keys fuel (Nat.gcd s k)
    | none => s

/-- Exact sign of `val x`, or `none` if the fuel (recursion depth) runs out.
A rigorous rational enclosure settles almost every case; otherwise
`x = A + B √k` and the sign follows from `sign A`, `sign B` and
`sign (A² - k B²)` (`rsign_split`). -/
def sign : ℕ → MQ → Option ℤ
  | 0, _ => none
  | fuel + 1, x =>
    let x := norm x
    if 0 < (iv x).1 then some 1
    else if (iv x).2 < 0 then some (-1)
    else match x.find? (fun p => 1 < p.1) with
      | none => some (qsign (ratPart x))
      | some p =>
        let k := refine (x.map (·.1)) p.1 p.1
        if k = 0 then none
        else
          match sign fuel (split k x).1, sign fuel (split k x).2 with
          | some sa, some sb =>
            if sb = 0 then some sa
            else if sa = 0 then some sb
            else if sa = sb then some sa
            else (sign fuel (sub (mul (split k x).1 (split k x).1)
                    (smul k (mul (split k x).2 (split k x).2)))).map (sa * ·)
          | _, _ => none

theorem sign_sound : ∀ (fuel : ℕ) (x : MQ) (s : ℤ), sign fuel x = some s → s = rsign (val x)
  | 0, x, s, h => by simp [sign] at h
  | fuel + 1, x, s, h => by
    have ih := sign_sound fuel
    have hv := val_norm x
    have hiv := iv_sound (norm x)
    unfold sign at h
    try simp only at h
    rw [← hv]
    by_cases h1 : 0 < (iv (norm x)).1
    · rw [if_pos h1] at h
      cases h
      exact (rsign_of_pos (lt_of_lt_of_le (by exact_mod_cast h1) hiv.1)).symm
    rw [if_neg h1] at h
    by_cases h2 : (iv (norm x)).2 < 0
    · rw [if_pos h2] at h
      cases h
      exact (rsign_of_neg (lt_of_le_of_lt hiv.2 (by exact_mod_cast h2))).symm
    rw [if_neg h2] at h
    split at h
    · rename_i hnone
      cases h
      rw [val_ratPart _ (fun p hp => by simpa using List.find?_eq_none.mp hnone p hp), qsign_eq]
    rename_i p _
    by_cases hk : refine ((norm x).map (·.1)) p.1 p.1 = 0
    · simp [hk] at h
    rw [if_neg hk] at h
    set k := refine ((norm x).map (·.1)) p.1 p.1
    have hk1 : (0 : ℝ) < k := by exact_mod_cast Nat.pos_of_ne_zero hk
    have hs := val_split k (norm x)
    have hsq := Real.sqrt_nonneg (k : ℝ)
    set A := (split k (norm x)).1
    set B := (split k (norm x)).2
    cases ha : sign fuel A with
    | none => simp [ha] at h
    | some sa =>
    cases hb : sign fuel B with
    | none => simp [ha, hb] at h
    | some sb =>
    simp only [ha, hb] at h
    have ea := ih A sa ha
    have eb := ih B sb hb
    by_cases c1 : sb = 0
    · rw [if_pos c1] at h
      cases h
      have hB : val B = 0 := rsign_eq_zero.mp (eb.symm.trans c1)
      rw [hs, hB, zero_mul, add_zero, ea]
    rw [if_neg c1] at h
    by_cases c2 : sa = 0
    · rw [if_pos c2] at h
      cases h
      have hA : val A = 0 := rsign_eq_zero.mp (ea.symm.trans c2)
      rw [hs, hA, zero_add, eb, rsign_mul_pos (Real.sqrt_pos.mpr hk1)]
    rw [if_neg c2] at h
    by_cases c3 : sa = sb
    · rw [if_pos c3] at h
      cases h
      rw [hs, ea]
      rcases rsign_cases (val A) with ⟨hA, e⟩ | ⟨hA, e⟩ | ⟨hA, e⟩
      · have hB : 0 < val B := by
          rcases rsign_cases (val B) with ⟨hB, _⟩ | ⟨_, f⟩ | ⟨_, f⟩
          · exact hB
          all_goals omega
        rw [e, rsign_of_pos (by nlinarith)]
      · have hB : val B < 0 := by
          rcases rsign_cases (val B) with ⟨_, f⟩ | ⟨hB, _⟩ | ⟨_, f⟩
          · omega
          · exact hB
          · omega
        rw [e, rsign_of_neg (by nlinarith)]
      · omega
    rw [if_neg c3] at h
    obtain ⟨sd, hd, rfl⟩ := Option.map_eq_some_iff.mp h
    have ed := ih _ sd hd
    rw [val_sub, val_smul, val_mul, val_mul] at ed
    rw [hs, ea, ed, rsign_split hk1
      (fun h0 => c2 (by rw [ea, h0, rsign_zero]))
      (fun h0 => c1 (by rw [eb, h0, rsign_zero]))
      (by rw [← ea, ← eb]; exact c3)]
    congr 2
    push_cast
    ring

end MQ

end VLink
