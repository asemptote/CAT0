import VLink.Kernel

/-!
# Completeness of the bead enumeration (docs/proofs.md §3)

Every bead the kernel accepts (`checkBead`) is found by its own enumeration
(`enumerateBeads`): `enumerateBeads_complete`. A π-bead is found when
π-beads are enumerated, as a bead with the same strip ends (its point on the
half circle is not recorded); any other bead is found exactly.

The search grows strips from each corner, keeping the interval of directions
at the start point `S` that cross every edge so far. The direction of an
accepted bead's geodesic lies strictly inside the interval of each edge it
crosses (`visible_crossing`), hence inside their intersections
(`intersect_within`), so the search never abandons the bead's strip.
-/

namespace VLink

/-! ### Determinants with a common first row -/

lemma rdet_swap23 (S X Y : Fin 3 → ℝ) : rdet S Y X = -rdet S X Y := by
  simp only [rdet_eq]; ring

/-- The Plücker relation for the 2-form `(x, y) ↦ det(S, x, y)`. -/
lemma rdet_plucker (S x y u v : Fin 3 → ℝ) :
    rdet S x y * rdet S u v = rdet S x u * rdet S y v - rdet S x v * rdet S y u := by
  simp only [rdet_eq]; ring

lemma rdet_of_cross_zero {S X : Fin 3 → ℝ} (h : crossProduct S X = 0) (Y : Fin 3 → ℝ) :
    rdet S X Y = 0 := by
  rw [cross_apply] at h
  have h0 := congrFun h 0
  have h1 := congrFun h 1
  have h2 := congrFun h 2
  simp at h0 h1 h2
  simp only [rdet_eq]
  linear_combination Y 0 * h0 + Y 1 * h1 + Y 2 * h2

lemma parallel_zero {G : Gram} {S X : V3} {p : Int} (h : parallel G S X = .ok p)
    {Y : Fin 3 → ℝ} (hd : rdet S.val X.val Y ≠ 0) : p = 0 := by
  rcases parallel_ok h with ⟨h0, -⟩ | ⟨-, hp⟩
  · exact absurd (rdet_of_cross_zero h0 Y) hd
  · exact hp

/-! ### Crossing an edge, and intervals of directions -/

/-- The geodesic leaving `S` in the direction of `w` meets the open arc `AB`
(`S` off the great circle of `AB`). -/
def Crossing (S w A B : Fin 3 → ℝ) : Prop :=
  rdet S A B ≠ 0 ∧ ∃ α β γ δ : ℝ, 0 < α ∧ 0 < β ∧ 0 < δ ∧ α • A + β • B = γ • S + δ • w

lemma Crossing.swap {S w A B : Fin 3 → ℝ} (h : Crossing S w A B) : Crossing S w B A := by
  obtain ⟨hd, α, β, γ, δ, hα, hβ, hδ, hc⟩ := h
  refine ⟨by rw [rdet_swap23]; exact neg_ne_zero.mpr hd, β, α, γ, δ, hβ, hα, hδ, ?_⟩
  rw [← hc, add_comm]

lemma Crossing.dets {S w A B : Fin 3 → ℝ} (h : Crossing S w A B) :
    0 < rdet S A w * rdet S A B ∧ 0 < rdet S w B * rdet S A B := by
  obtain ⟨hd, α, β, γ, δ, hα, hβ, hδ, hc⟩ := h
  have e1 : δ * rdet S A w = β * rdet S A B := by
    have := congrArg (fun v => rdet S A v) hc
    simp only [rdet_eq, Pi.add_apply, Pi.smul_apply, smul_eq_mul] at this ⊢
    linear_combination -this
  have e2 : δ * rdet S w B = α * rdet S A B := by
    have := congrArg (fun v => rdet S v B) hc
    simp only [rdet_eq, Pi.add_apply, Pi.smul_apply, smul_eq_mul] at this ⊢
    linear_combination -this
  have hD2 : 0 < rdet S A B * rdet S A B := mul_self_pos.mpr hd
  have f1 : δ * (rdet S A w * rdet S A B) = β * (rdet S A B * rdet S A B) := by
    rw [← mul_assoc, e1]; ring
  have f2 : δ * (rdet S w B * rdet S A B) = α * (rdet S A B * rdet S A B) := by
    rw [← mul_assoc, e2]; ring
  constructor
  · by_contra hn
    push Not at hn
    nlinarith [mul_pos hβ hD2, mul_nonneg hδ.le (neg_nonneg.mpr hn)]
  · by_contra hn
    push Not at hn
    nlinarith [mul_pos hα hD2, mul_nonneg hδ.le (neg_nonneg.mpr hn)]

/-- `w` points strictly inside the interval `(L, R)` of directions at `S`
(an interval narrower than `π`). -/
def Within (S w : Fin 3 → ℝ) (iv : Vec × Vec) : Prop :=
  0 < rdet S iv.1.val iv.2.val ∧ 0 < rdet S iv.1.val w ∧ 0 < rdet S w iv.2.val

/-- `X` lies in the closed interval `I` of directions at `S`. -/
def InI (S : Fin 3 → ℝ) (I : Vec × Vec) (X : Fin 3 → ℝ) : Prop :=
  0 ≤ rdet S I.1.val X ∧ 0 ≤ rdet S X I.2.val

/-- The visible-interval lemma: a geodesic from `S` crossing the open arc
`AB` leaves `S` strictly inside the interval `visible` computes. -/
lemma visible_crossing {S A B : Vec} {w : Fin 3 → ℝ} {o : Option (Vec × Vec)}
    (hc : Crossing S.val w A.val B.val) (h : visible S A B = .ok o) :
    ∃ iv, o = some iv ∧ Within S.val w iv := by
  obtain ⟨d1, d2⟩ := hc.dets
  have hD := hc.1
  unfold visible at h
  split at h
  · cases h
  rename_i s hs
  cases h
  rw [sgn_ok hs, val_det3]
  rcases lt_or_gt_of_ne hD with hneg | hpos
  · rw [rsign_of_neg hneg]
    refine ⟨(B, A), by simp, ?_⟩
    have x1 : rdet S.val A.val w < 0 := by
      by_contra hn; push Not at hn
      nlinarith [mul_nonneg hn (le_of_lt (neg_pos.mpr hneg))]
    have x2 : rdet S.val w B.val < 0 := by
      by_contra hn; push Not at hn
      nlinarith [mul_nonneg hn (le_of_lt (neg_pos.mpr hneg))]
    refine ⟨?_, ?_, ?_⟩
    · show 0 < rdet S.val B.val A.val
      rw [rdet_swap23]; linarith
    · show 0 < rdet S.val B.val w
      rw [rdet_swap23]; linarith
    · show 0 < rdet S.val w A.val
      rw [rdet_swap23]; linarith
  · rw [rsign_of_pos hpos]
    refine ⟨(A, B), by simp, hpos, ?_, ?_⟩
    · show 0 < rdet S.val A.val w
      by_contra hn; push Not at hn
      nlinarith [mul_nonneg (neg_nonneg.mpr hn) hpos.le]
    · show 0 < rdet S.val w B.val
      by_contra hn; push Not at hn
      nlinarith [mul_nonneg (neg_nonneg.mpr hn) hpos.le]

lemma inside_strict {S L R X : Vec} {c : Bool} (hL : 0 < rdet S.val L.val X.val)
    (hR : 0 < rdet S.val X.val R.val) (h : inside S L R X true = .ok c) : c = true := by
  unfold inside at h
  split at h
  · rename_i a b ha hb
    cases h
    rw [sgn_ok ha, sgn_ok hb, val_det3, val_det3, rsign_of_pos hL, rsign_of_pos hR]
    decide
  · cases h
  · cases h

lemma inside_closed {S L R X : Vec} {c : Bool} (h : inside S L R X false = .ok c) :
    (c = true ↔ InI S.val (L, R) X.val) := by
  unfold inside at h
  split at h
  · rename_i a b ha hb
    cases h
    rw [sgn_ok ha, sgn_ok hb, val_det3, val_det3]
    simp only [Bool.false_eq_true, ↓reduceIte, decide_eq_true_eq, InI]
    rw [ge_iff_le, ge_iff_le, rsign_nonneg, rsign_nonneg]
  · cases h
  · cases h

lemma pickEnd_spec {S : Vec} {I1 I2 : Vec × Vec} {X1 X2 : Vec} {o : Option Vec}
    (h : pickEnd S I1 I2 X1 X2 = .ok o) :
    (o = some X2 ∧ InI S.val I1 X2.val) ∨ (o = some X1 ∧ InI S.val I2 X1.val) ∨
      (o = none ∧ ¬ InI S.val I1 X2.val ∧ ¬ InI S.val I2 X1.val) := by
  unfold pickEnd at h
  split at h
  · cases h
  · rename_i h1
    cases h
    exact Or.inl ⟨rfl, (inside_closed h1).mp rfl⟩
  · rename_i h1
    have n1 : ¬ InI S.val I1 X2.val := fun hi => by simpa using (inside_closed h1).mpr hi
    split at h
    · cases h
    · rename_i h2
      cases h
      exact Or.inr (Or.inl ⟨rfl, (inside_closed h2).mp rfl⟩)
    · rename_i h2
      cases h
      exact Or.inr (Or.inr ⟨rfl, n1, fun hi => by simpa using (inside_closed h2).mpr hi⟩)

/-- The intersection of two intervals containing `w` contains `w`. -/
lemma intersect_within {S : Vec} {w : Fin 3 → ℝ} {I1 I2 : Vec × Vec} {o : Option (Vec × Vec)}
    (h1 : Within S.val w I1) (h2 : Within S.val w I2) (h : intersect S I1 I2 = .ok o) :
    ∃ new, o = some new ∧ Within S.val w new := by
  obtain ⟨L1, R1⟩ := I1
  obtain ⟨L2, R2⟩ := I2
  obtain ⟨a1, b1, c1⟩ := h1
  obtain ⟨a2, b2, c2⟩ := h2
  simp only at a1 b1 c1 a2 b2 c2
  set s := S.val
  have P1 := rdet_plucker s L2.val R1.val L1.val w
  have P2 := rdet_plucker s L1.val R2.val L2.val w
  have P3 := rdet_plucker s L1.val R2.val w R1.val
  have P4 := rdet_plucker s L2.val R1.val w R2.val
  rw [rdet_swap23 s L1.val L2.val, rdet_swap23 s w R1.val, rdet_swap23 s L1.val R1.val] at P1
  rw [rdet_swap23 s L2.val L1.val, rdet_swap23 s w R2.val, rdet_swap23 s L2.val R2.val] at P2
  rw [rdet_swap23 s w R2.val] at P3
  rw [rdet_swap23 s w R1.val] at P4
  have F1 : 0 ≤ rdet s L1.val L2.val → 0 < rdet s L2.val R1.val := by
    intro hh; by_contra hn; push Not at hn
    nlinarith [mul_nonneg hh c1.le, mul_pos b2 a1]
  have F2 : 0 ≤ rdet s L2.val L1.val → 0 < rdet s L1.val R2.val := by
    intro hh; by_contra hn; push Not at hn
    nlinarith [mul_nonneg hh c2.le, mul_pos b1 a2]
  have F3 : 0 ≤ rdet s R2.val R1.val → 0 < rdet s L1.val R2.val := by
    intro hh; by_contra hn; push Not at hn
    nlinarith [mul_nonneg b1.le hh, mul_pos a1 c2]
  have F4 : 0 ≤ rdet s R1.val R2.val → 0 < rdet s L2.val R1.val := by
    intro hh; by_contra hn; push Not at hn
    nlinarith [mul_nonneg b2.le hh, mul_pos a2 c1]
  have sw12 := rdet_swap23 s L1.val L2.val
  have swR := rdet_swap23 s R1.val R2.val
  unfold intersect at h
  simp only at h
  rcases hLs : pickEnd S (L1, R1) (L2, R2) L1 L2 with e | oL
  · rw [hLs] at h; cases h
  rcases hRs : pickEnd S (L1, R1) (L2, R2) R1 R2 with e | oR
  · rw [hLs, hRs] at h; simp at h
  -- the left end is `L2` or `L1`
  obtain ⟨L, rfl, hLc⟩ : ∃ L, oL = some L ∧
      ((L = L2 ∧ 0 ≤ rdet s L1.val L2.val) ∨ (L = L1 ∧ 0 ≤ rdet s L2.val L1.val)) := by
    rcases pickEnd_spec hLs with ⟨e, hi⟩ | ⟨e, hi⟩ | ⟨e, n1, n2⟩
    · exact ⟨L2, e, Or.inl ⟨rfl, hi.1⟩⟩
    · exact ⟨L1, e, Or.inr ⟨rfl, hi.1⟩⟩
    · exfalso
      by_cases hh : 0 ≤ rdet s L1.val L2.val
      · exact n1 ⟨hh, (F1 hh).le⟩
      · exact n2 ⟨by linarith, (F2 (by linarith)).le⟩
  obtain ⟨R, rfl, hRc⟩ : ∃ R, oR = some R ∧
      ((R = R2 ∧ 0 ≤ rdet s R2.val R1.val) ∨ (R = R1 ∧ 0 ≤ rdet s R1.val R2.val)) := by
    rcases pickEnd_spec hRs with ⟨e, hi⟩ | ⟨e, hi⟩ | ⟨e, n1, n2⟩
    · exact ⟨R2, e, Or.inl ⟨rfl, hi.2⟩⟩
    · exact ⟨R1, e, Or.inr ⟨rfl, hi.2⟩⟩
    · exfalso
      by_cases hh : 0 ≤ rdet s R2.val R1.val
      · exact n1 ⟨(F3 hh).le, hh⟩
      · exact n2 ⟨(F4 (by linarith)).le, by linarith⟩
  have hLR : 0 < rdet s L.val R.val := by
    rcases hLc with ⟨rfl, hh⟩ | ⟨rfl, hh⟩ <;> rcases hRc with ⟨rfl, -⟩ | ⟨rfl, -⟩
    · exact a2
    · exact F1 hh
    · exact F2 hh
    · exact a1
  have hLw : 0 < rdet s L.val w := by rcases hLc with ⟨rfl, -⟩ | ⟨rfl, -⟩ <;> assumption
  have hwR : 0 < rdet s w R.val := by rcases hRc with ⟨rfl, -⟩ | ⟨rfl, -⟩ <;> assumption
  rw [hLs, hRs] at h
  simp only at h
  split at h
  · rename_i sv hsv
    cases h
    rw [sgn_ok hsv, val_det3, rsign_of_pos hLR]
    exact ⟨(L, R), by simp, hLR, hLw, hwR⟩
  · cases h

/-! ### The enumeration only adds beads -/

/-- Every bead of `o` is in `o'`. -/
def Grows (o o' : Array Bead) : Prop := ∀ x ∈ o, x ∈ o'

lemma Grows.refl (o : Array Bead) : Grows o o := fun _ h => h

lemma Grows.trans {o1 o2 o3 : Array Bead} (h1 : Grows o1 o2) (h2 : Grows o2 o3) : Grows o1 o3 :=
  fun x hx => h2 x (h1 x hx)

lemma Grows.push (o : Array Bead) (x : Bead) : Grows o (o.push x) :=
  fun _ hy => Array.mem_push.mpr (Or.inl hy)

lemma recordBead_grows {G : Gram} {pb : Bool} {ts2 : Array Nat} {pos2 : Array (Array Vec)}
    {steps2 : Array Nat} {start r : Nat} {S E : Vec} {iv : Vec × Vec} {par : Int}
    {out out' : Array Bead}
    (h : recordBead G pb ts2 pos2 steps2 start r S E iv par out = .ok out') : Grows out out' := by
  unfold recordBead at h
  split_ifs at h
  · split at h
    · cases h; exact Grows.push _ _
    · cases h
  · cases h; exact Grows.refl _
  · split at h
    · cases h
    · cases h; exact Grows.refl _
    · split at h
      · cases h; exact Grows.push _ _
      · cases h
  · cases h; exact Grows.refl _

lemma growBranch_grows {G : Gram} {rec : Nat → Vec × Vec → Array Bead → KM (Array Bead)}
    {S E : Vec} {P2 : Array Vec} {r dn : Nat} {iv : Vec × Vec} {out out' : Array Bead}
    (h : growBranch G rec S E P2 r dn iv out = .ok out')
    (hrec : ∀ dn new o o', rec dn new o = .ok o' → Grows o o') : Grows out out' := by
  unfold growBranch at h
  split at h
  · cases h
  · split_ifs at h
    · cases h; exact Grows.refl _
    · split at h
      · cases h
      · cases h; exact Grows.refl _
      · split at h
        · cases h
        · cases h; exact Grows.refl _
        · exact hrec _ _ _ _ h

lemma grow_grows {L : Link} {G : Gram} {pb : Bool} :
    ∀ (fuel : Nat) (ts : Array Nat) (pos : Array (Array Vec)) (steps : Array Nat) (start : Nat)
      (S : Vec) (d : Nat) (iv : Vec × Vec) (out out' : Array Bead),
      grow L G pb fuel ts pos steps start S d iv out = .ok out' → Grows out out' := by
  intro fuel
  induction fuel with
  | zero => intro _ _ _ _ _ _ _ _ _ h; simp [grow] at h
  | succ fuel ih =>
    intro ts pos steps start S d iv out out' h
    unfold grow at h
    dsimp only at h
    split at h
    · cases h
    · split at h
      · cases h
      · rename_i out1 hrec
        have g1 := recordBead_grows hrec
        split_ifs at h
        · cases h; exact g1
        · split at h
          · cases h
          · rename_i out2 hb
            exact g1.trans
              ((growBranch_grows hb (fun _ _ _ _ h' => ih _ _ _ _ _ _ _ _ _ h')).trans
                (growBranch_grows h (fun _ _ _ _ h' => ih _ _ _ _ _ _ _ _ _ h')))

lemma enumStart_grows {L : Link} {pb : Bool} {t d0 : Nat} {out out' : Array Bead}
    (h : enumStart L pb t d0 out = .ok out') : Grows out out' := by
  unfold enumStart at h
  dsimp only at h
  split at h
  · cases h
  rename_i ba _
  split at h
  · cases h
  rename_i bb _
  have g0 : Grows out ((out.push ba).push bb) := (Grows.push out _).trans (Grows.push _ _)
  split at h
  · split_ifs at h
    · cases h; exact g0
    · split at h
      · cases h
      · cases h; exact g0
      · exact g0.trans (grow_grows _ _ _ _ _ _ _ _ _ _ h)
  · cases h
  · cases h

lemma enumList_ok {L : Link} {pb : Bool} :
    ∀ (l : List (Nat × Nat)) (out res : Array Bead), enumList L pb l out = .ok res →
      Grows out res ∧ ∀ p ∈ l, ∃ o o', enumStart L pb p.1 p.2 o = .ok o' ∧ Grows o' res := by
  intro l
  induction l with
  | nil => intro out res h; cases h; exact ⟨Grows.refl _, by simp⟩
  | cons p l ih =>
    intro out res h
    unfold enumList at h
    split at h
    · cases h
    · rename_i out' hs
      obtain ⟨g, hall⟩ := ih _ _ h
      refine ⟨(enumStart_grows hs).trans g, ?_⟩
      intro q hq
      rcases List.mem_cons.mp hq with rfl | hq
      · exact ⟨out, out', hs, g⟩
      · exact hall q hq

/-! ### Arrays -/

lemma getElem!_push_lt {α : Type} [Inhabited α] (a : Array α) (x : α) {i : Nat} (hi : i < a.size) :
    (a.push x)[i]! = a[i]! := by
  rw [getElem!_pos (a.push x) i (by simp; omega), getElem!_pos a i hi, Array.getElem_push_lt]

lemma getElem!_push_eq {α : Type} [Inhabited α] (a : Array α) (x : α) {i : Nat} (hi : a.size = i) :
    (a.push x)[i]! = x := by
  subst hi
  rw [getElem!_pos (a.push x) a.size (by simp), Array.getElem_push_eq]

lemma back!_of_size {α : Type} [Inhabited α] {a : Array α} {k : Nat} (h : a.size = k + 1) :
    a.back! = a[k]! := by
  unfold Array.back!; rw [h]; rfl

lemma array_ext! {α : Type} [Inhabited α] {a b : Array α} (h1 : a.size = b.size)
    (h2 : ∀ i < a.size, a[i]! = b[i]!) : a = b := by
  apply Array.ext h1
  intro i hi hi'
  have := h2 i hi
  rwa [getElem!_pos a i hi, getElem!_pos b i hi'] at this

/-! ### Facts about an accepted bead -/

lemma Bead.make_pi {ts : Array Nat} {start stop : Nat} {G : Gram} {pos : Array (Array Vec)}
    {steps : Array Nat} {Q : Option Vec} {b : Bead} {par : Int}
    (h : Bead.make ts start stop G pos steps Q = .ok b)
    (hp : parallel G (pos[0]!)[start]! pos.back![stop]! = .ok par) : b.pi = decide (par < 0) := by
  unfold Bead.make at h
  dsimp only at h
  rw [hp] at h
  simp only [Except.ok.injEq] at h
  rw [← h]

lemma develop_zero {L : Link} {ts steps : Array Nat} {pos : Array (Array V3)}
    (h : develop L ts steps = .ok pos) : pos[0]! = E3 := by
  unfold develop at h
  split_ifs at h
  split at h
  · cases h; simp
  · cases h

/-- What `checkBead` decided about the bead's endpoints. -/
lemma checkBead_par {L : Link} {spec : BeadSpec} {b : Bead} (h : checkBead L spec = .ok b) :
    b.G = L.G[spec.ts[0]!]! ∧ spec.stop < 3 ∧ ∃ par, parallel b.G b.S b.E = .ok par ∧ par ≤ 0 ∧
      b.pi = decide (par < 0) ∧
      (par = 0 → Bead.make spec.ts spec.start spec.stop b.G b.pos spec.steps none = .ok b ∧
        ¬ (spec.ts.size = 1 ∧ spec.start = spec.stop)) := by
  unfold checkBead at h
  obtain ⟨_, hm, h⟩ := bind_ok h
  obtain ⟨pos, -, h⟩ := bind_ok h
  obtain ⟨_, -, h⟩ := bind_ok h
  obtain ⟨par, hpar, h⟩ := bind_ok h
  obtain ⟨_, hgd, h⟩ := bind_ok h
  obtain ⟨Q, hQ, h⟩ := bind_ok h
  obtain ⟨b', hb, h⟩ := bind_ok h
  obtain ⟨_, -, h⟩ := bind_ok h
  cases h
  obtain ⟨-, hG, hposb, -, hs, he, hS, hE, -⟩ := Bead.make_ok hb
  have hm' := guardK_ok hm
  simp only [Bool.not_eq_true', Bool.or_eq_false_iff, decide_eq_false_iff_not, not_lt,
    Array.any_eq_false, not_le] at hm'
  have hg' := guardK_ok hgd
  refine ⟨hG, by omega, par, ?_, by simpa using hg', Bead.make_pi hb hpar, ?_⟩
  · rw [hG, hS, hE]; exact hpar
  · intro h0
    subst h0
    simp only [lt_irrefl, ↓reduceIte] at hQ
    split_ifs at hQ with c2
    cases hQ
    refine ⟨by rw [hG, hposb]; exact hb, ?_⟩
    rintro ⟨h1, h2⟩
    apply c2
    simp [h1, h2]

/-- The direction in which an accepted bead's geodesic leaves `S`: toward
`E`, or for a π-bead along the half circle through its point (with the sign
that makes the crossing coefficients positive). -/
noncomputable def beadDir (b : Bead) : Fin 3 → ℝ :=
  if b.pi then
    (if 0 < bform b.G.toMat b.outDir.val b.outDir.val * bform b.G.toMat b.S.val b.S.val -
        bform b.G.toMat b.outDir.val b.S.val ^ 2 then (1 : ℝ) else -1) • b.outDir.val
  else b.E.val

/-- An accepted bead's geodesic crosses every interior edge of its strip. -/
lemma bead_crossing {L : Link} {spec : BeadSpec} {b : Bead} (h : checkBead L spec = .ok b)
    (hsym : ∀ u v, bform b.G.toMat u v = bform b.G.toMat v u) :
    ∀ k, k + 1 < b.ts.size → Crossing b.S.val (beadDir b)
      ((b.pos[k]!)[(others3 b.steps[k]!).1]!).val ((b.pos[k]!)[(others3 b.steps[k]!).2]!).val := by
  obtain ⟨-, -, -, -, -, -, -, -, -, -, -, hout, hcr⟩ := checkBead_ok h
  intro k hk
  have hc := hcr k (List.mem_range.mpr (by omega))
  cases hpi : b.pi
  · rw [hpi, hout hpi] at hc
    obtain ⟨α, β, γ, δ, hα, hβ, hγ, hδ, hcomb, hdS⟩ := crosses_ok_minor hc
    refine ⟨by rw [rdet_cycle]; exact hdS, α, β, γ, δ, hα, hβ, hδ, ?_⟩
    simp only [beadDir, hpi, Bool.false_eq_true, ↓reduceIte]
    exact hcomb
  · rw [hpi] at hc
    obtain ⟨α, β, γ, δ, hα, hβ, hcomb, hK, hdS⟩ := crosses_ok_half hsym hc
    refine ⟨by rw [rdet_cycle]; exact hdS, ?_⟩
    simp only [beadDir, hpi, ↓reduceIte]
    split_ifs with hK0
    · refine ⟨α, β, γ, δ, hα, hβ, ?_, by rw [one_smul]; exact hcomb⟩
      by_contra hn; push Not at hn
      nlinarith [mul_nonneg (neg_nonneg.mpr hn) hK0.le]
    · refine ⟨α, β, γ, -δ, hα, hβ, ?_, ?_⟩
      · by_contra hn; push Not at hn
        push Not at hK0
        nlinarith [mul_nonneg (neg_nonneg.mpr hn) (neg_nonneg.mpr hK0)]
      · rw [hcomb, smul_smul]; congr 1; ring_nf

/-- A bead found by the enumeration for the accepted bead `b`: the same bead,
or for a π-bead one with the same strip and ends. -/
def BeadMatch (b b' : Bead) : Prop :=
  b'.pi = b.pi ∧ b'.ts = b.ts ∧ b'.start = b.start ∧ b'.stop = b.stop ∧ (b.pi = false → b' = b)

/-- The search state after the first `k + 1` triangles of `b`'s strip, about to
cross edge `k`: the developed strip so far, and an interval of directions
strictly containing the bead's. -/
structure Lvl (b : Bead) (k : Nat) (tsE : Array Nat) (posE : Array (Array Vec))
    (stepsE : Array Nat) (d : Nat) (iv : Vec × Vec) : Prop where
  hk : k + 1 < b.ts.size
  ts_size : tsE.size = k + 1
  ts_eq : ∀ i < k + 1, tsE[i]! = b.ts[i]!
  pos_size : posE.size = k + 1
  pos_eq : ∀ i < k + 1, posE[i]! = b.pos[i]!
  steps_size : stepsE.size = k
  steps_eq : ∀ i < k, stepsE[i]! = b.steps[i]!
  hd : d = b.steps[k]!
  within : Within b.S.val (beadDir b) iv

lemma growBranch_hit {G : Gram} {rec : Nat → Vec × Vec → Array Bead → KM (Array Bead)}
    {S E : Vec} {P2 : Array Vec} {r dn : Nat} {iv : Vec × Vec} {out out' : Array Bead}
    {w : Fin 3 → ℝ} (hc : Crossing S.val w E.val (P2[3 - r - dn]!).val) (hiv : Within S.val w iv)
    (h : growBranch G rec S E P2 r dn iv out = .ok out') :
    ∃ new, Within S.val w new ∧ rec dn new out = .ok out' := by
  unfold growBranch at h
  split at h
  · cases h
  rename_i p hp
  have hp0 : p = 0 := parallel_zero hp (Y := E.val) (by rw [rdet_swap23]; exact neg_ne_zero.mpr hc.1)
  subst hp0
  rw [if_neg (by simp)] at h
  split at h
  · cases h
  · rename_i hv
    obtain ⟨_, e, -⟩ := visible_crossing hc hv
    cases e
  · rename_i iv2 hv
    obtain ⟨_, e, hw2⟩ := visible_crossing hc hv
    cases e
    split at h
    · cases h
    · rename_i hi
      obtain ⟨_, e, -⟩ := intersect_within hiv hw2 hi
      cases e
    · rename_i new hi
      obtain ⟨_, e, hw⟩ := intersect_within hiv hw2 hi
      cases e
      exact ⟨new, hw, h⟩

/-- A strip's end vertex is not an end of the edge just crossed. -/
lemma end_off_edge {S w A B E : Fin 3 → ℝ} (hc : Crossing S w A B) (hE : E = A ∨ E = B)
    (hw : w = E ∨ crossProduct S E = 0) : False := by
  have dA : rdet S A A = 0 := rdet_self_right S A
  have dB : rdet S B B = 0 := rdet_self_right S B
  obtain ⟨d1, d2⟩ := hc.dets
  rcases hE with rfl | rfl <;> rcases hw with rfl | hx
  · rw [dA, zero_mul] at d1; exact lt_irrefl _ d1
  · exact hc.1 (rdet_of_cross_zero hx B)
  · rw [dB, zero_mul] at d2; exact lt_irrefl _ d2
  · apply hc.1
    rw [← neg_neg (rdet S A E), ← rdet_swap23, rdet_of_cross_zero hx A, neg_zero]

/-- The search follows an accepted bead's strip from any state on it to the
end, and records the bead. -/
theorem grow_complete {geo : Geometry} {iv : Nat} {L : Link} (hL : Link.make geo iv = .ok L)
    {spec : BeadSpec} {b : Bead} (h : checkBead L spec = .ok b) {piBeads : Bool}
    (hpi : piBeads = true ∨ b.pi = false) :
    ∀ (fuel k : Nat) (tsE : Array Nat) (posE : Array (Array Vec)) (stepsE : Array Nat) (d : Nat)
      (ivl : Vec × Vec) (out out' : Array Bead), Lvl b k tsE posE stepsE d ivl →
      grow L b.G piBeads fuel tsE posE stepsE b.start b.S d ivl out = .ok out' →
      ∃ b' ∈ out', BeadMatch b b' := by
  obtain ⟨hts, hst, hs, he, hn, -, -, hdev, hnb, hS, hE, hout, -⟩ := checkBead_ok h
  obtain ⟨hG, hstop, parB, hparB, hparB0, hpiB, hmake⟩ := checkBead_par h
  rw [← hts, ← hst] at hdev
  obtain ⟨hpsz, hsteps, hdstep⟩ := develop_step hdev
  have hn' : ∀ i < b.ts.size, b.ts[i]! < L.n := by
    intro i hi; apply hn; rw [← hts, getElem!_pos b.ts i hi]; simp
  have hsym : ∀ u v, bform b.G.toMat u v = bform b.G.toMat v u := by
    have h0 : spec.ts[0]! < L.n := by rw [← hts]; exact hn' 0 (by omega)
    rw [hG]
    exact bform_symm (mkTrans_t2 (Link.make_trans hL _ h0 0 (by norm_num))).2.1
  have hcross := bead_crossing h hsym
  have hstop' : b.stop < 3 := by rw [he]; exact hstop
  intro fuel
  induction fuel with
  | zero => intro k _ _ _ _ _ _ _ _ hg; simp [grow] at hg
  | succ fuel ih =>
    intro k tsE posE stepsE d ivl out out' hlv hg
    obtain ⟨hk, htsz, htse, hpsz', hpse, hssz, hsse, hd, hwin⟩ := hlv
    subst hd
    have hks : k < b.steps.size := by omega
    obtain ⟨hPk, hdk, hadj⟩ := hdstep k hks
    have htk := hn' k (by omega)
    obtain ⟨hg', hperm, -⟩ := mkTrans_tr (Link.make_trans hL _ htk _ hdk)
    obtain ⟨iA, iB, -⟩ := stepWith_idx (L.tr b.ts[k]! b.steps[k]!) b.pos[k]! hdk hperm
    obtain ⟨ha3, hb3, hab, had, hbd⟩ := others3_spec hdk
    obtain ⟨gab3, gab⟩ := isPerm3_spec hperm ha3 hb3 hab
    obtain ⟨-, gad⟩ := isPerm3_spec hperm ha3 hdk had
    obtain ⟨gb3, gbd⟩ := isPerm3_spec hperm hb3 hdk hbd
    obtain ⟨gd3, -⟩ := isPerm3_spec hperm hdk ha3 (Ne.symm had)
    have hbt : tsE.back! = b.ts[k]! := by rw [back!_of_size htsz]; exact htse k (by omega)
    have hbp : posE.back! = b.pos[k]! := by rw [back!_of_size hpsz']; exact hpse k (by omega)
    have hP2 : (step L posE.back! tsE.back! b.steps[k]!).2 = b.pos[k + 1]! := by
      rw [hbt, hbp, hPk]; rfl
    have hr : ((L.glue[tsE.back!]!)[b.steps[k]!]!)[b.steps[k]!]! =
        (L.tr b.ts[k]! b.steps[k]!).g[b.steps[k]!]! := by rw [hbt, hg']
    have hadj' : (L.adj[tsE.back!]!)[b.steps[k]!]! = b.ts[k + 1]! := by rw [hbt, hadj]
    unfold grow at hg
    dsimp only at hg
    rw [hP2, hr, hadj'] at hg
    -- the state after crossing edge `k`
    have nts : ∀ i < k + 2, (tsE.push b.ts[k + 1]!)[i]! = b.ts[i]! := by
      intro i hi
      by_cases hi' : i < k + 1
      · rw [getElem!_push_lt _ _ (by omega)]; exact htse i hi'
      · rw [getElem!_push_eq _ _ (by omega), show i = k + 1 by omega]
    have npos : ∀ i < k + 2, (posE.push b.pos[k + 1]!)[i]! = b.pos[i]! := by
      intro i hi
      by_cases hi' : i < k + 1
      · rw [getElem!_push_lt _ _ (by omega)]; exact hpse i hi'
      · rw [getElem!_push_eq _ _ (by omega), show i = k + 1 by omega]
    have nsteps : ∀ i < k + 1, (stepsE.push b.steps[k]!)[i]! = b.steps[i]! := by
      intro i hi
      by_cases hi' : i < k
      · rw [getElem!_push_lt _ _ (by omega)]; exact hsse i hi'
      · rw [getElem!_push_eq _ _ (by omega), show i = k by omega]
    have hcross_k := hcross k hk
    set g := (L.tr b.ts[k]! b.steps[k]!).g with hgdef
    set r := g[b.steps[k]!]! with hrdef
    split at hg
    · cases hg
    rename_i par hpar
    split at hg
    · cases hg
    rename_i out1 hrec
    -- everything after recording only adds
    have hafter : (∃ b' ∈ out1, BeadMatch b b') → ∃ b' ∈ out', BeadMatch b b' := by
      rintro ⟨b', hb', hm⟩
      refine ⟨b', ?_, hm⟩
      split_ifs at hg
      · cases hg; exact hb'
      · split at hg
        · cases hg
        · rename_i out2 hb2
          exact growBranch_grows hg (fun _ _ _ _ h' => grow_grows _ _ _ _ _ _ _ _ _ _ h') _
            (growBranch_grows hb2 (fun _ _ _ _ h' => grow_grows _ _ _ _ _ _ _ _ _ _ h') _ hb')
    by_cases hlast : k + 2 = b.ts.size
    · -- the last edge: the strip ends at the new vertex, and the bead is recorded
      have hbpos : b.pos.back! = b.pos[k + 1]! := back!_of_size (by omega)
      have hstopr : b.stop = r := by
        by_contra hne
        have hcases : b.stop = g[(others3 b.steps[k]!).1]! ∨
            b.stop = g[(others3 b.steps[k]!).2]! := by
          omega
        have hEe : b.E.val = ((b.pos[k]!)[(others3 b.steps[k]!).1]!).val ∨
            b.E.val = ((b.pos[k]!)[(others3 b.steps[k]!).2]!).val := by
          rw [hE, hbpos, hPk]
          rcases hcases with hc | hc
          · left; rw [hc, iA]
          · right; rw [hc, iB]
        refine end_off_edge hcross_k hEe ?_
        cases hpb : b.pi
        · left; simp [beadDir, hpb]
        · right
          rcases parallel_ok hparB with ⟨hx, -⟩ | ⟨-, hx⟩
          · exact hx
          · rw [hpiB, hx] at hpb; simp at hpb
      have hEP : (b.pos[k + 1]!)[r]! = b.E := by rw [hE, hbpos, hstopr]
      rw [hEP] at hpar hrec
      have hpp : par = parB := by rw [hparB] at hpar; cases hpar; rfl
      subst hpp
      have hts2 : tsE.push b.ts[k + 1]! = b.ts :=
        array_ext! (by simp [htsz]; omega) (fun i hi => nts i (by simp [htsz] at hi; omega))
      have hpos2 : posE.push b.pos[k + 1]! = b.pos :=
        array_ext! (by simp [hpsz']; omega) (fun i hi => npos i (by simp [hpsz'] at hi; omega))
      have hsteps2 : stepsE.push b.steps[k]! = b.steps :=
        array_ext! (by simp [hssz]; omega) (fun i hi => nsteps i (by simp [hssz] at hi; omega))
      rw [hts2, hpos2, hsteps2, ← hstopr] at hrec
      apply hafter
      unfold recordBead at hrec
      by_cases hneg : par < 0
      · rw [if_pos hneg] at hrec
        have hpiT : b.pi = true := by rw [hpiB]; simpa using hneg
        have hpb : piBeads = true := by
          rcases hpi with h1 | h1
          · exact h1
          · rw [hpiT] at h1; cases h1
        rw [if_pos hpb] at hrec
        split at hrec
        · rename_i b'' hb''
          cases hrec
          obtain ⟨t1, -, -, -, s1, s2, -⟩ := Bead.make_ok hb''
          have hp'' := Bead.make_pi hb'' (par := par) (by rw [← hS, ← hE]; exact hparB)
          refine ⟨b'', Array.mem_push.mpr (Or.inr rfl), ?_, t1, s1, s2, ?_⟩
          · rw [hp'', hpiT]; simpa using hneg
          · intro hf; rw [hpiT] at hf; cases hf
        · cases hrec
      · rw [if_neg hneg] at hrec
        have hz : par = 0 := by omega
        rw [if_pos hz] at hrec
        have hpiF : b.pi = false := by rw [hpiB]; simp [hz]
        have hdir : beadDir b = b.E.val := by simp [beadDir, hpiF]
        rw [hdir] at hwin
        split at hrec
        · cases hrec
        · rename_i hin
          have := inside_strict hwin.2.1 hwin.2.2 hin
          cases this
        · split at hrec
          · rename_i b'' hb''
            cases hrec
            have hmk := (hmake hz).1
            rw [← hts, ← hs, ← he, ← hst] at hmk
            rw [hmk] at hb''
            cases hb''
            exact ⟨b, Array.mem_push.mpr (Or.inr rfl), rfl, rfl, rfl, rfl, fun _ => rfl⟩
          · cases hrec
    · -- an interior edge follows: the search continues into the right branch
      have hk2 : k + 2 < b.ts.size := by omega
      obtain ⟨-, hdk1, -⟩ := hdstep (k + 1) (by omega)
      have hnb' : b.steps[k + 1]! ≠ r := by
        have := noBacktrack_spec hnb (k + 1) (by omega) (by rw [← hst]; omega)
        simp only [Nat.add_sub_cancel] at this
        have e : r = ((L.glue[b.ts[k]!]!)[b.steps[k]!]!)[b.steps[k]!]! := by
          rw [hrdef]; exact congrArg (fun a : Array Nat => a[b.steps[k]!]!) hg'
        rw [e]
        simpa [hts, hst] using this
      set dn := b.steps[k + 1]! with hdn
      have hy3 : 3 - r - dn < 3 := by omega
      have hc1 := hcross (k + 1) hk2
      have hcE : Crossing b.S.val (beadDir b) ((b.pos[k + 1]!)[r]!).val
          ((b.pos[k + 1]!)[3 - r - dn]!).val := by
        rcases others3_pair hdk1 gd3 hy3 (Ne.symm hnb') (by omega) (by omega) with hp | hp
        · rw [hp] at hc1; exact hc1
        · rw [hp] at hc1; exact hc1.swap
      have hp0 : par = 0 := parallel_zero hpar hcE.1
      rw [if_neg (by simp [hp0])] at hg
      have mk : ∀ new, Within b.S.val (beadDir b) new →
          Lvl b (k + 1) (tsE.push b.ts[k + 1]!) (posE.push b.pos[k + 1]!)
            (stepsE.push b.steps[k]!) dn new := fun new hw =>
        ⟨hk2, by simp [htsz], nts, by simp [hpsz'], npos, by simp [hssz], nsteps, rfl, hw⟩
      split at hg
      · cases hg
      rename_i out2 hb1
      rcases (show dn = (r + 1) % 3 ∨ dn = (r + 2) % 3 by omega) with h1 | h1
      · rw [← h1] at hb1
        obtain ⟨new, hw, hrc⟩ := growBranch_hit hcE hwin hb1
        obtain ⟨b', hb', hm⟩ := ih (k + 1) _ _ _ _ _ _ _ (mk new hw) hrc
        exact ⟨b', growBranch_grows hg (fun _ _ _ _ h' => grow_grows _ _ _ _ _ _ _ _ _ _ h') _ hb', hm⟩
      · rw [← h1] at hg
        obtain ⟨new, hw, hrc⟩ := growBranch_hit hcE hwin hg
        exact ih (k + 1) _ _ _ _ _ _ _ (mk new hw) hrc

lemma E3_val (i : Nat) (hi : i < 3) (j : Fin 3) :
    (E3[i]!).val j = if (j : Nat) = i then 1 else 0 := by
  interval_cases i <;> fin_cases j <;> simp [E3, V3.val, mqOne, MQ.ofRat]

lemma E3_cross {i j : Nat} (hi : i < 3) (hj : j < 3)
    (h : crossProduct (E3[i]!).val (E3[j]!).val = 0) : i = j := by
  rw [cross_apply] at h
  have h0 := congrFun h 0
  have h1 := congrFun h 1
  have h2 := congrFun h 2
  simp only [E3_val i hi, E3_val j hj] at h0 h1 h2
  interval_cases i <;> interval_cases j <;> simp_all

lemma E3_ne {i : Nat} (hi : i < 3) : (E3[i]!).val ≠ 0 := by
  intro h
  have := congrFun h ⟨i, hi⟩
  rw [E3_val i hi] at this
  simp at this

lemma rdet_self_mid (S X : Fin 3 → ℝ) : rdet S S X = 0 := by simp only [rdet_eq]; ring

/-- The search from the corner where an accepted bead starts records it. -/
theorem enumStart_complete {geo : Geometry} {iv : Nat} {L : Link} (hL : Link.make geo iv = .ok L)
    {spec : BeadSpec} {b : Bead} (h : checkBead L spec = .ok b) {piBeads : Bool}
    (hpi : piBeads = true ∨ b.pi = false) {out out' : Array Bead}
    (hen : enumStart L piBeads b.ts[0]! b.start out = .ok out') : ∃ b' ∈ out', BeadMatch b b' := by
  obtain ⟨hts, hst, hs, he, hn, hsz, hstart, hdev, -, hS, hE, -, -⟩ := checkBead_ok h
  obtain ⟨hG, hstop, parB, hparB, hparB0, -, hmake⟩ := checkBead_par h
  rw [← hts, ← hst] at hdev
  obtain ⟨hpsz, hsteps, hdstep⟩ := develop_step hdev
  have hp0 := develop_zero hdev
  have hstart' : b.start < 3 := by rw [hs]; exact hstart
  have hstop' : b.stop < 3 := by rw [he]; exact hstop
  have hsz' : b.ts.size ≠ 0 := by rw [hts]; exact hsz
  have ht0 : b.ts[0]! < L.n := by
    apply hn; rw [← hts, getElem!_pos b.ts 0 (by omega)]; simp
  have hG' : L.G[b.ts[0]!]! = b.G := by rw [hG, hts]
  have hS' : E3[b.start]! = b.S := by rw [hS, hp0]
  have hsym : ∀ u v, bform b.G.toMat u v = bform b.G.toMat v u := by
    rw [← hG']
    exact bform_symm (mkTrans_t2 (Link.make_trans hL _ ht0 0 (by norm_num))).2.1
  unfold enumStart at hen
  dsimp only at hen
  split at hen
  · cases hen
  rename_i ba hba
  split at hen
  · cases hen
  rename_i bb hbb
  have hrest : Grows ((out.push ba).push bb) out' := by
    split at hen
    · split_ifs at hen
      · cases hen; exact Grows.refl _
      · split at hen
        · cases hen
        · cases hen; exact Grows.refl _
        · exact grow_grows _ _ _ _ _ _ _ _ _ _ hen
    · cases hen
    · cases hen
  by_cases h1 : b.ts.size = 1
  · -- an edge of the triangle
    have hpos1 : b.pos.size = 1 := by omega
    have hEe : b.E = E3[b.stop]! := by rw [hE, back!_of_size hpos1, hp0]
    have hpar0 : parB = 0 := by
      rcases parallel_ok hparB with ⟨hx, hr⟩ | ⟨-, hr⟩
      · exfalso
        rw [← hS', hEe] at hx
        have hij := E3_cross hstart' hstop' hx
        have hpos := link_gram_posdef hL ht0 _ (E3_ne hstart')
        rw [hr, ← hG', ← hS', hEe, ← hij, rsign_of_pos hpos] at hparB0
        omega
      · exact hr
    obtain ⟨hmk, hnd⟩ := hmake hpar0
    have hne : b.start ≠ b.stop := fun e =>
      hnd ⟨by rw [← hts]; exact h1, by rw [← hs, ← he]; exact e⟩
    have hbts : b.ts = #[b.ts[0]!] := array_ext! (by simp [h1]) (fun i hi => by
      have : i = 0 := by omega
      subst this; simp)
    have hbsteps : b.steps = #[] := array_ext! (by show b.steps.size = 0; omega) (fun i hi => by omega)
    have hbpos : b.pos = #[E3] := array_ext! (by simp [hpos1]) (fun i hi => by
      have : i = 0 := by omega
      subst this; simp [hp0])
    rw [← hts, ← hst, ← hs, ← he, hbsteps, hbpos, ← hG'] at hmk
    have hmk' : Bead.make #[b.ts[0]!] b.start b.stop L.G[b.ts[0]!]! #[E3] #[] none = .ok b := by
      rw [← hbts]; exact hmk
    have hmem : b ∈ (out.push ba).push bb := by
      rcases (show b.stop = (b.start + 1) % 3 ∨ b.stop = (b.start + 2) % 3 by omega) with hq | hq
      · rw [hq] at hmk'
        rw [hmk'] at hba
        cases hba
        exact Array.mem_push.mpr (Or.inl (Array.mem_push.mpr (Or.inr rfl)))
      · rw [hq] at hmk'
        rw [hmk'] at hbb
        cases hbb
        exact Array.mem_push.mpr (Or.inr rfl)
    exact ⟨b, hrest b hmem, rfl, rfl, rfl, rfl, fun _ => rfl⟩
  · -- a longer strip: the search follows it
    have h2 : 1 < b.ts.size := by omega
    have hc0 := bead_crossing h hsym 0 h2
    have hd0lt : b.steps[0]! < 3 := (hdstep 0 (by omega)).2.1
    have hd0 : b.steps[0]! = b.start := by
      by_contra hne
      have hS0 : b.S.val = ((b.pos[0]!)[b.start]!).val := by rw [hS]
      rcases others3_mem hd0lt hstart' (Ne.symm hne) with hm | hm
      · apply hc0.1
        rw [hS0, hm]; exact rdet_self_mid _ _
      · apply hc0.1
        rw [hS0, hm]; exact rdet_self_left _ _
    rw [hd0, hp0] at hc0
    have hcab : Crossing b.S.val (beadDir b) (E3[(b.start + 1) % 3]!).val
        (E3[(b.start + 2) % 3]!).val := by
      rcases others3_pair hstart' (p := (b.start + 1) % 3) (q := (b.start + 2) % 3)
        (by omega) (by omega) (by omega) (by omega) (by omega) with hp | hp
      · rw [hp] at hc0; exact hc0
      · rw [hp] at hc0; exact hc0.swap
    rw [← hS'] at hcab
    split at hen
    · rename_i pa pb hpa hpb
      have pa0 : pa = 0 := parallel_zero hpa hcab.1
      have pb0 : pb = 0 := parallel_zero hpb hcab.swap.1
      rw [if_neg (by simp [pa0, pb0])] at hen
      split at hen
      · cases hen
      · rename_i hv
        obtain ⟨_, e, -⟩ := visible_crossing hcab hv
        cases e
      · rename_i iv0 hv
        obtain ⟨_, e, hw⟩ := visible_crossing hcab hv
        cases e
        rw [hS'] at hw
        rw [hG', hS'] at hen
        refine grow_complete hL h hpi _ 0 _ _ _ _ _ _ _
          ⟨by omega, rfl, ?_, rfl, ?_, rfl, ?_, hd0.symm, hw⟩ hen
        · intro i hi
          have : i = 0 := by omega
          subst this; simp
        · intro i hi
          have : i = 0 := by omega
          subst this; simp [hp0]
        · intro i hi; omega
    · cases hen
    · cases hen

/-- **Completeness of the bead enumeration.** Every bead the kernel accepts
is found by `enumerateBeads` (π-beads when they are enumerated). -/
theorem enumerateBeads_complete {geo : Geometry} {iv : Nat} {L : Link} {piBeads : Bool}
    {beads : Array Bead} (hL : Link.make geo iv = .ok L)
    (hen : enumerateBeads L piBeads = .ok beads)
    {spec : BeadSpec} {b : Bead} (h : checkBead L spec = .ok b)
    (hpi : piBeads = true ∨ b.pi = false) :
    ∃ i < beads.size, BeadMatch b beads[i]! := by
  obtain ⟨hts, -, hs, -, hn, hsz, hstart, -⟩ := checkBead_ok h
  obtain ⟨-, hall⟩ := enumList_ok _ _ _ hen
  have hmem : (b.ts[0]!, b.start) ∈ corners L := by
    have ht0 : b.ts[0]! < L.n := by
      have hsz' : 0 < b.ts.size := by rw [hts]; omega
      apply hn; rw [← hts, getElem!_pos b.ts 0 hsz']; simp
    simp only [corners, List.mem_flatMap, List.mem_range, List.mem_map, Prod.mk.injEq]
    exact ⟨b.ts[0]!, ht0, b.start, by rw [hs]; exact hstart, rfl, rfl⟩
  obtain ⟨o, o', hst, hg⟩ := hall _ hmem
  obtain ⟨b', hb', hm⟩ := enumStart_complete hL h hpi hst
  obtain ⟨i, hi, e⟩ := Array.mem_iff_getElem.mp (hg b' hb')
  exact ⟨i, hi, by rw [getElem!_pos beads i hi, e]; exact hm⟩

end VLink
