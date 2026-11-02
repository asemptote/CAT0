import VLink.Cert
import VLink.Complete

/-!
# What an accepted certificate proves (layer L4)

`LocallyCAT0 isosig sqlens`: the piecewise Euclidean metric on the closed
triangulated 3-manifold `isosig` with squared edge lengths `sqlens` is locally
CAT(0). Mathlib has no piecewise Euclidean complexes, so the predicate is
opaque, and the theorems from the literature that characterise it
(`docs/proofs.md` §0) enter as the axioms below, each stated in terms of what
the kernel computes:

* `not_cat0_of_degenerate`, `not_cat0_of_short_edge`: Gromov's link condition
  (Bridson–Haefliger Thm II.5.2) and circles (BH Example II.3.17): an edge's link is a circle of
  length its dihedral sum;
* `not_cat0_of_necklace`: (G), (C) and local geodesics at cone points (BH I.5.7, I.7.39):
  a short necklace is a closed local geodesic of length `< 2π` in a vertex link;
* `cat0_of_no_short_necklace`: (G), (C), Bowditch's shrinking (EMM Lemma 4.5)
  and the gallery reduction (Elder–McCammond 2002; EMM Lemmas 4.7–4.8): if a
  vertex link is not CAT(1), it has a short necklace (of beads the kernel
  accepts, with junctions `≥ π`);
* `cat0_of_no_short_necklace_nopi`: the same with π-beads left out, relying on
  EMM 2003 Lemma 4.9 as well.

So the axioms say which kernel objects the cited theorems are about; both
directions use the same `ShortNecklace`. Everything else is proved: the exact
arithmetic and angle sums, the meaning of each geometric predicate,
development, bead and junction checks, the completeness of the kernel's bead
enumeration (`VLink/Complete.lean`: every accepted bead is enumerated), and
the shortest-path certificates. Hence

* `checkFalse_sound`: an accepted False certificate ⇒ `¬ LocallyCAT0`;
* `checkTrue_sound`: an accepted True certificate ⇒ `LocallyCAT0`.

`#print axioms checkTrue_sound` lists the cited theorems it uses, next to
Lean's standard axioms.
-/

namespace VLink

open Real

/-- The metric is locally CAT(0) (opaque: characterised by the axioms below). -/
opaque LocallyCAT0 (isosig : String) (sqlens : Array ℚ) : Prop

/-- The total dihedral angle around edge `e`. -/
noncomputable def edgeSum (geo : Geometry) (e : Nat) : ℝ :=
  ((geo.edgeAngles e).toList.map fun p => angleOf p.1 p.2).sum

/-- The length of a bead (`π` for a π-bead). -/
noncomputable def beadAngle (b : Bead) : ℝ := angleOf (beadLength b).1 (beadLength b).2

/-- A necklace may go from bead `x` to bead `y`: `y` starts where `x` ends and
the junction passes (both angles `≥ π`; counted as passing when a π-bead is
involved, which only adds necklaces). -/
def Passes (L : Link) (beads : Array Bead) (x y : Nat) : Prop :=
  x < beads.size ∧ ∃ row, junctionRow L beads x = .ok row ∧ (y, true) ∈ row

/-- A necklace of length `< 2π` in the link of vertex `iv`: beads accepted by
the kernel (geodesics crossing their strips' interior edges; π-beads only if
`piBeads`), meeting one after the other, with both junction angles (sums of
corner angles) `≥ π`. -/
def ShortNecklace (geo : Geometry) (iv : Nat) (piBeads : Bool) : Prop :=
  ∃ L : Link, Link.make geo iv = .ok L ∧ ∃ beads : List Bead, beads ≠ [] ∧
    (∀ b ∈ beads, (∃ s, checkBead L s = .ok b) ∧ (piBeads = true ∨ b.pi = false)) ∧
    (beads.map beadAngle).sum < 2 * π ∧
    ∀ k < beads.length, endV L beads[k]! = startV L beads[(k + 1) % beads.length]! ∧
      ∀ side ∈ [1, 2], ∃ r, junctionFan L beads[k]! beads[(k + 1) % beads.length]! side = .ok r ∧
        π ≤ (r.2.map fun p => angleOf p.1 p.2).sum

/-- No necklace of the enumerated beads is shorter than `2π` (what a True
certificate's trees show). -/
def NoShortNecklace (geo : Geometry) (iv : Nat) (piBeads : Bool) : Prop :=
  ∃ L : Link, Link.make geo iv = .ok L ∧ ∃ beads : Array Bead, enumerateBeads L piBeads = .ok beads ∧
    (∀ x < beads.size, ∃ row, junctionRow L beads x = .ok row) ∧
    ∀ l : List ℕ, l ≠ [] → (∀ i ∈ l, i < beads.size) → (l : Cycle ℕ).Chain (Passes L beads) →
      2 * π ≤ (l.map fun i => beadAngle beads[i]!).sum

/-! ### The cited theorems -/

/-- A degenerate tetrahedron carries no Euclidean metric: by convention the
metric is not locally CAT(0) (the certificate reason "invalid"). -/
axiom not_cat0_of_degenerate {isosig : String} {sqlens : Array ℚ} {T : Tri} {geo : Geometry} :
    Tri.ofIsoSig isosig = .ok T → Geometry.make T sqlens = .ok geo →
    ∀ t < T.n, geo.tetValid t = false → ¬ LocallyCAT0 isosig sqlens

/-- (G), (C): the link of an edge is a circle of length its dihedral sum,
CAT(1) iff that is `≥ 2π`. -/
axiom not_cat0_of_short_edge {isosig : String} {sqlens : Array ℚ} {T : Tri} {geo : Geometry} :
    Tri.ofIsoSig isosig = .ok T → Geometry.make T sqlens = .ok geo →
    (∀ t < T.n, geo.tetValid t = true) →
    ∀ e < T.edges.size, edgeSum geo e < 2 * π → ¬ LocallyCAT0 isosig sqlens

/-- (G), (C), (L): a short necklace is a closed local geodesic of length
`< 2π` in a vertex link, so that link is not CAT(1). -/
axiom not_cat0_of_necklace {isosig : String} {sqlens : Array ℚ} {T : Tri} {geo : Geometry} :
    Tri.ofIsoSig isosig = .ok T → Geometry.make T sqlens = .ok geo →
    (∀ t < T.n, geo.tetValid t = true) →
    ∀ iv < T.corners.size, ShortNecklace geo iv true → ¬ LocallyCAT0 isosig sqlens

/-- (G), (C), (B), (E): valid tetrahedra, every edge `≥ 2π` and no short
necklace at any vertex make the metric locally CAT(0). -/
axiom cat0_of_no_short_necklace {isosig : String} {sqlens : Array ℚ} {T : Tri} {geo : Geometry} :
    Tri.ofIsoSig isosig = .ok T → Geometry.make T sqlens = .ok geo →
    (∀ t < T.n, geo.tetValid t = true) → (∀ e < T.edges.size, 2 * π ≤ edgeSum geo e) →
    (∀ iv < T.corners.size, ¬ ShortNecklace geo iv true) → LocallyCAT0 isosig sqlens

/-- The same with beads of length exactly `π` left out (EMM 2003, Lemma 4.9). -/
axiom cat0_of_no_short_necklace_nopi {isosig : String} {sqlens : Array ℚ} {T : Tri}
    {geo : Geometry} :
    Tri.ofIsoSig isosig = .ok T → Geometry.make T sqlens = .ok geo →
    (∀ t < T.n, geo.tetValid t = true) → (∀ e < T.edges.size, 2 * π ≤ edgeSum geo e) →
    (∀ iv < T.corners.size, ¬ ShortNecklace geo iv false) → LocallyCAT0 isosig sqlens

/-! ### Soundness -/

lemma atLeast_rep {st : AngleSum} {Θ : ℝ} (hr : st.Rep Θ) (k : Nat) :
    atLeast st k = true ↔ (k : ℝ) * π ≤ Θ := by
  rw [← AngleSum.atLeast_iff hr k]
  simp [atLeast]

/-! ### From the enumerated beads to all accepted beads -/

lemma forall2_mem_left {α β : Type} {R : α → β → Prop} {l : List α} {l' : List β}
    (h : List.Forall₂ R l l') {a : α} (ha : a ∈ l) : ∃ b ∈ l', R a b := by
  induction h with
  | nil => simp at ha
  | cons hab _ ih =>
    rcases List.mem_cons.mp ha with rfl | ha
    · exact ⟨_, List.mem_cons_self, hab⟩
    · obtain ⟨b, hb, hR⟩ := ih ha
      exact ⟨b, List.mem_cons_of_mem _ hb, hR⟩

lemma beadAngle_match {b b' : Bead} (h : BeadMatch b b') : beadAngle b' = beadAngle b := by
  obtain ⟨hpi, -, -, -, heq⟩ := h
  cases hb : b.pi
  · rw [heq hb]
  · simp [beadAngle, beadLength, hpi, hb]

lemma startV_match {L : Link} {b b' : Bead} (h : BeadMatch b b') : startV L b' = startV L b := by
  obtain ⟨-, hts, hs, -, -⟩ := h
  simp [startV, hts, hs]

lemma endV_match {L : Link} {b b' : Bead} (h : BeadMatch b b') : endV L b' = endV L b := by
  obtain ⟨-, hts, -, he, -⟩ := h
  simp [endV, hts, he]

lemma junctionOk_of_fans {L : Link} {bin bout : Bead}
    (h : ∀ side ∈ [1, 2], ∃ r, junctionFan L bin bout side = .ok r ∧
      π ≤ (r.2.map fun p => angleOf p.1 p.2).sum) : junctionOk L bin bout = .ok true := by
  obtain ⟨r1, hr1, h1⟩ := h 1 (by simp)
  obtain ⟨r2, hr2, h2⟩ := h 2 (by simp)
  have a1 : atLeast r1.1 1 = true := (atLeast_rep (junctionFan_ok hr1) 1).mpr (by simpa using h1)
  have a2 : atLeast r2.1 1 = true := (atLeast_rep (junctionFan_ok hr2) 1).mpr (by simpa using h2)
  have s1 : junctionSide L bin bout 1 = .ok r1.1 := by simp [junctionSide, hr1]
  have s2 : junctionSide L bin bout 2 = .ok r2.1 := by simp [junctionSide, hr2]
  unfold junctionOk
  rw [s1, s2]
  show Except.ok (atLeast r1.1 1 && atLeast r2.1 1) = _
  rw [a1, a2]
  rfl

/-- If the enumerated beads have no short necklace, no beads the kernel
accepts do: every accepted bead is enumerated (`enumerateBeads_complete`). -/
theorem not_short_of_enum {geo : Geometry} {iv : Nat} {piBeads : Bool}
    (h : NoShortNecklace geo iv piBeads) : ¬ ShortNecklace geo iv piBeads := by
  rintro ⟨L', hL', bs, hne, hchk, hlen, hjunc⟩
  obtain ⟨L, hL, beads, hen, hrows, hall⟩ := h
  have hLL : L' = L := by rw [hL] at hL'; exact (Except.ok.inj hL').symm
  rw [hLL] at hchk hjunc
  set n := bs.length with hn
  have hn0 : n ≠ 0 := by simpa [hn] using hne
  -- each bead of the necklace is an enumerated bead
  have hfind : ∀ k, ∃ i, k < n → i < beads.size ∧ BeadMatch bs[k]! beads[i]! := by
    intro k
    by_cases hk : k < n
    · have hmem : bs[k]! ∈ bs := by rw [getElem!_pos bs k hk]; exact List.getElem_mem _
      obtain ⟨⟨s, hs⟩, hp⟩ := hchk _ hmem
      obtain ⟨i, hi, hm⟩ := enumerateBeads_complete hL hen hs hp
      exact ⟨i, fun _ => ⟨hi, hm⟩⟩
    · exact ⟨0, fun h => absurd h hk⟩
  choose idx hidx using hfind
  -- consecutive enumerated beads pass
  have hpass : ∀ k < n, Passes L beads (idx k) (idx ((k + 1) % n)) := by
    intro k hk
    have hk1 : (k + 1) % n < n := Nat.mod_lt _ (by omega)
    obtain ⟨hix, hmx⟩ := hidx k hk
    obtain ⟨hiy, hmy⟩ := hidx _ hk1
    obtain ⟨hmeet, hfans⟩ := hjunc k hk
    obtain ⟨row, hrow⟩ := hrows _ hix
    refine ⟨hix, row, hrow, ?_⟩
    unfold junctionRow at hrow
    have hy : idx ((k + 1) % n) ∈
        (List.range beads.size).filter (fun y => startV L beads[y]! == endV L beads[idx k]!) := by
      simp only [List.mem_filter, List.mem_range, beq_iff_eq]
      refine ⟨hiy, ?_⟩
      rw [startV_match hmy, endV_match hmx]
      exact hmeet.symm
    obtain ⟨e, he, hfe⟩ := forall2_mem_left (mapK_ok hrow) hy
    have hentry : junctionEntry L beads (idx k) (idx ((k + 1) % n)) =
        .ok (idx ((k + 1) % n), true) := by
      unfold junctionEntry
      split_ifs with hpi
      · rfl
      · simp only [Bool.or_eq_true, not_or, Bool.not_eq_true] at hpi
        rw [hmx.2.2.2.2 (by rw [← hmx.1]; exact hpi.1),
          hmy.2.2.2.2 (by rw [← hmy.1]; exact hpi.2), junctionOk_of_fans hfans]
    rw [hentry] at hfe
    cases hfe
    exact he
  obtain ⟨n', hn'⟩ : ∃ n', n = n' + 1 := ⟨n - 1, by omega⟩
  have hchain : (((List.range n).map idx : List ℕ) : Cycle ℕ).Chain (Passes L beads) := by
    rw [← Cycle.map_coe, Cycle.chain_map, hn', Cycle.chain_range_succ]
    refine ⟨?_, fun m hm => ?_⟩
    · have := hpass n' (by omega)
      rwa [show (n' + 1) % n = 0 by rw [hn']; exact Nat.mod_self _] at this
    · have := hpass m (by omega)
      rwa [show (m + 1) % n = m + 1 from Nat.mod_eq_of_lt (by omega)] at this
  have hge := hall ((List.range n).map idx) (by simp [hn0]) (by
    intro i hi
    obtain ⟨k, hk, rfl⟩ := List.mem_map.mp hi
    exact (hidx k (List.mem_range.mp hk)).1) hchain
  have hsum : (((List.range n).map idx).map fun i => beadAngle beads[i]!).sum =
      (bs.map beadAngle).sum := by
    rw [List.map_map]
    congr 1
    apply List.ext_getElem (by simp [hn])
    intro i h1 h2
    have hi : i < n := by simpa using h1
    simp only [List.getElem_map, List.getElem_range, Function.comp_apply]
    rw [beadAngle_match (hidx i hi).2, getElem!_pos bs i hi]
  rw [hsum] at hge
  push_cast at hge
  linarith


/-- An accepted False certificate proves the metric is not locally CAT(0). -/
theorem checkFalse_sound {c : FalseCert} {m : String} (h : checkFalse c = .ok m) :
    ¬ LocallyCAT0 c.isosig c.sqlens := by
  unfold checkFalse at h
  split at h
  · cases h
  rename_i T hT
  split at h
  · cases h
  rename_i geo hg
  split at h
  · -- a degenerate tetrahedron
    rename_i t _
    split_ifs at h with hbad
    simp only [Bool.or_eq_true, decide_eq_true_eq, not_or, Bool.not_eq_true] at hbad
    exact not_cat0_of_degenerate hT hg t (by omega) hbad.2
  · -- a short edge
    rename_i e _
    split at h
    · cases h
    rename_i hv
    split_ifs at h with he
    split at h
    · cases h
    rename_i acc hacc
    split_ifs at h with hat
    unfold Geometry.edgeAngle at hacc
    have hrep := angleSum_ok hacc
    refine not_cat0_of_short_edge hT hg (tetsValidK_ok hv) e (by omega) ?_
    by_contra hle
    push Not at hle
    exact hat ((atLeast_rep hrep 2).mpr (by push_cast; exact hle))
  · -- a short necklace
    rename_i iv necklace _
    split at h
    · cases h
    rename_i hv
    split_ifs at h with hiv
    split at h
    · cases h
    rename_i L hL
    split at h
    · cases h
    rename_i beads hbeads
    split_ifs at h with hemp
    split at h
    · cases h
    rename_i tot htot
    split_ifs at h with hat
    split at h
    · cases h
    rename_i hneck
    refine not_cat0_of_necklace hT hg (tetsValidK_ok hv) iv (by omega) ⟨L, hL, beads, ?_, ?_, ?_, ?_⟩
    · intro h0; simp [h0] at hemp
    · intro b hb
      obtain ⟨s, -, hsb⟩ := forall2_mem_right (mapK_ok hbeads) hb
      exact ⟨⟨s, hsb⟩, Or.inl rfl⟩
    · unfold totalLength at htot
      have hrep := angleSum_ok htot
      by_contra hle
      push Not at hle
      apply hat
      refine (atLeast_rep hrep 2).mpr ?_
      have hsum : ((beads.toArray.map beadLength).toList.map fun p => angleOf p.1 p.2).sum =
          (beads.map beadAngle).sum := by
        simp [Function.comp_def] <;> rfl
      rw [hsum]
      push_cast
      linarith
    · intro k hk
      have hk' := forallK_ok hneck k (List.mem_range.mpr hk)
      split_ifs at hk' with hmeet
      refine ⟨not_not.mp hmeet, ?_⟩
      intro side hside
      have hj := forallK_ok hk' side hside
      split at hj
      · rename_i acc hacc
        split_ifs at hj with hat1
        unfold junctionSide at hacc
        split at hacc
        · rename_i r hr
          cases hacc
          refine ⟨r, hr, ?_⟩
          have := (atLeast_rep (junctionFan_ok hr) 1).mp hat1
          simpa using this
        · cases hacc
      · cases hj

/-- An accepted True certificate proves the metric is locally CAT(0). -/
theorem checkTrue_sound {c : TrueCert} {m : String} (h : checkTrue c = .ok m) :
    LocallyCAT0 c.isosig c.sqlens := by
  unfold checkTrue at h
  split_ifs at h with h49 hrep
  split at h
  · cases h
  rename_i T hT
  split at h
  · cases h
  rename_i geo hg
  split at h
  · cases h
  rename_i hv
  split at h
  · cases h
  rename_i hedges
  dsimp only at h
  split_ifs at h with hcover
  split at h
  · cases h
  rename_i hverts
  have hvalid := tetsValidK_ok hv
  have hedge : ∀ e < T.edges.size, 2 * π ≤ edgeSum geo e := by
    intro e he
    have := forallK_ok hedges e (List.mem_range.mpr he)
    split at this
    · rename_i acc hacc
      split_ifs at this with hat
      unfold Geometry.edgeAngle at hacc
      have := (atLeast_rep (angleSum_ok hacc) 2).mp hat
      push_cast at this
      exact this
    · cases this
  have hvert : ∀ iv < T.corners.size, NoShortNecklace geo iv c.piBeads := by
    intro iv hiv
    have hmem : iv ∈ (c.vertices.qsort (fun a b => a.1 < b.1)).toList.map (·.1) := by
      rw [not_not.mp hcover]; exact List.mem_range.mpr hiv
    obtain ⟨v, hvm, hv1⟩ := List.mem_map.mp hmem
    have hok := forallK_ok hverts v hvm
    unfold vertexOkK at hok
    rw [hv1] at hok
    split at hok
    · cases hok
    rename_i L hL
    split at hok
    · cases hok
    rename_i beads hb
    split_ifs at hok with hn
    obtain ⟨Jt, hrows, hcyc⟩ := checkVertexTrees_ok hok
    refine ⟨L, hL, beads, hb, fun x hx => ⟨_, hrows x hx⟩, ?_⟩
    intro l hl hln hc
    apply hcyc l hl hln
    refine Cycle.Chain.imp ?_ hc
    rintro x y ⟨hx, row, hrow, hy⟩
    rw [hrows x hx] at hrow
    cases hrow
    exact hy
  have hno : ∀ iv < T.corners.size, ¬ ShortNecklace geo iv c.piBeads :=
    fun iv hiv => not_short_of_enum (hvert iv hiv)
  cases hpi : c.piBeads
  · exact cat0_of_no_short_necklace_nopi hT hg hvalid hedge (by simpa [hpi] using hno)
  · exact cat0_of_no_short_necklace hT hg hvalid hedge (by simpa [hpi] using hno)

end VLink
