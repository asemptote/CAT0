import Mathlib

/-!
# Shortest-path certificates rule out short necklaces (layer L4, pure)

The graph argument behind True certificates (`docs/true_certificates.md`).
Beads are `0, …, n-1` with lengths `len i ≥ 0`; `pass x y` says a necklace
may go from bead `x` to bead `y` (they meet and the junction passes). A
necklace is a cyclic sequence of beads, each passing to the next: a `Cycle`
whose `Chain` is `pass`.

For each bead `s` the certificate gives a set of tree nodes `T s` (beads
`≥ s`) with values `d s x` such that

* `s ∈ T s` and `d s s ≤ len s`;
* for a node `x` and a bead `y > s` with `pass x y`: either
  `d s x + len y ≥ 2π`, or `y` is a node with `d s y ≤ d s x + len y`;
* for a node `x` with `pass x s`: `d s x ≥ 2π`.

Then every necklace has length `≥ 2π` (`no_short_necklace`): walk around it
from its lowest bead `s`; by induction each prefix ends at a node whose value
is at most the prefix length, until the walk returns to `s`.
-/

namespace VLink

open Real

variable {n : ℕ} {len : ℕ → ℝ} {pass : ℕ → ℕ → Prop} {T : ℕ → ℕ → Prop} {d : ℕ → ℕ → ℝ}

/-- The walk from a node `x` (reached with prefix length `P ≥ d s x`) along
`rest`, then back to `s`, has length `≥ 2π`. -/
lemma walk_bound (hlen : ∀ i, 0 ≤ len i) (s : ℕ)
    (h1 : ∀ x, T s x → ∀ y, pass x y → s < y → y < n →
      2 * π ≤ d s x + len y ∨ (T s y ∧ d s y ≤ d s x + len y))
    (h2 : ∀ x, T s x → pass x s → 2 * π ≤ d s x) :
    ∀ (rest : List ℕ) (x : ℕ) (P : ℝ), T s x → d s x ≤ P →
      List.IsChain pass (x :: (rest ++ [s])) → (∀ y ∈ rest, s ≤ y ∧ y < n) →
      2 * π ≤ P + (rest.map len).sum := by
  intro rest
  induction rest with
  | nil =>
    intro x P hT hd hch _
    simp only [List.nil_append, List.isChain_pair] at hch
    simpa using (h2 x hT hch).trans hd
  | cons y r ih =>
    intro x P hT hd hch hge
    simp only [List.cons_append, List.isChain_cons_cons] at hch
    obtain ⟨hxy, hch⟩ := hch
    have hsum : 0 ≤ (r.map len).sum := List.sum_nonneg (by simp [hlen])
    obtain ⟨hsy, hyn⟩ := hge y List.mem_cons_self
    rcases hsy.lt_or_eq with hlt | heq
    · rcases h1 x hT y hxy hlt hyn with hfar | ⟨hTy, hdy⟩
      · simp only [List.map_cons, List.sum_cons]
        linarith [hlen y]
      · have := ih y (P + len y) hTy (by linarith) hch
          (fun z hz => hge z (List.mem_cons_of_mem _ hz))
        simp only [List.map_cons, List.sum_cons]
        linarith
    · subst heq
      have := h2 x hT hxy
      simp only [List.map_cons, List.sum_cons]
      linarith [hlen s]

/-- No necklace is shorter than `2π`. -/
theorem no_short_necklace (hlen : ∀ i, 0 ≤ len i)
    (h0 : ∀ s < n, T s s ∧ d s s ≤ len s)
    (h1 : ∀ s < n, ∀ x, T s x → ∀ y, pass x y → s < y → y < n →
      2 * π ≤ d s x + len y ∨ (T s y ∧ d s y ≤ d s x + len y))
    (h2 : ∀ s < n, ∀ x, T s x → pass x s → 2 * π ≤ d s x)
    (l : List ℕ) (hl : l ≠ []) (hn : ∀ i ∈ l, i < n) (hc : (l : Cycle ℕ).Chain pass) :
    2 * π ≤ (l.map len).sum := by
  -- the lowest bead s of the necklace
  obtain ⟨s, hs, hmin⟩ : ∃ s ∈ l, ∀ y ∈ l, s ≤ y := by
    obtain ⟨s, hs, hmin⟩ := l.toFinset.exists_min_image id
      (List.toFinset_nonempty_iff l |>.mpr hl)
    exact ⟨s, List.mem_toFinset.mp hs, fun y hy => hmin y (List.mem_toFinset.mpr hy)⟩
  -- rotate the necklace to start at s
  obtain ⟨a, b, rfl⟩ := List.append_of_mem hs
  have hrot : (a ++ s :: b) ~r (s :: (b ++ a)) := by
    have := List.isRotated_append (l := a) (l' := s :: b)
    simpa using this
  have hc' : ((s :: (b ++ a) : List ℕ) : Cycle ℕ).Chain pass := by
    rwa [← Cycle.coe_eq_coe.mpr hrot]
  rw [Cycle.chain_coe_cons] at hc'
  have hsum : ((a ++ s :: b).map len).sum = len s + ((b ++ a).map len).sum := by
    rw [(hrot.perm.map len).sum_eq]; simp
  rw [hsum]
  have hsn : s < n := hn s hs
  obtain ⟨hT0, hd0⟩ := h0 s hsn
  have hall : ∀ y ∈ b ++ a, s ≤ y ∧ y < n := by
    intro y hy
    have hy' : y ∈ a ++ s :: b := by
      rcases List.mem_append.mp hy with h | h <;> simp [h]
    exact ⟨hmin y hy', hn y hy'⟩
  have := walk_bound hlen s (h1 s hsn) (h2 s hsn) (b ++ a) s (len s) hT0 hd0 hc' hall
  linarith

end VLink
