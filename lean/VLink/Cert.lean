import VLink.Kernel

/-!
# Reading certificates (JSON, as written by `vlink/check.py` and `vlink/certify.py`)
-/

namespace VLink

open Lean (Json)

def FALSE_FORMAT : String := "vlink-certificate/2"
def TRUE_FORMAT : String := "vlink-certificate-true/2"

/-- `"p/q"` or `"p"`. -/
def parseRat (s : String) : Except String ℚ := do
  match s.splitOn "/" with
  | [p] => match p.trimAscii.toString.toInt? with
    | some p => return p
    | none => throw s!"bad rational {s}"
  | [p, q] => match p.trimAscii.toString.toInt?, q.trimAscii.toString.toNat? with
    | some p, some q => if q = 0 then throw "zero denominator" else return (p : ℚ) / q
    | _, _ => throw s!"bad rational {s}"
  | _ => throw s!"bad rational {s}"

def objList (j : Json) : Except String (List (String × Json)) := do
  return (← j.getObj?).toList

/-- A field element `{"s": "p/q", ...}` (keys: square-free radicands). -/
def parseMQ (j : Json) : Except String MQ := do
  let mut out : MQ := []
  for (k, v) in ← objList j do
    let some s := k.toNat? | throw s!"bad radicand {k}"
    let c : ℚ ← match v with
      | .str t => parseRat t
      | _ => do let i ← v.getInt?; pure (i : ℚ)
    out := out ++ [(s, c)]
  return out

def parseNats (j : Json) : Except String (Array Nat) := do
  (← j.getArr?).mapM (·.getNat?)

def parseSqlens (j : Json) : Except String (Array ℚ) := do
  (← j.getArr?).mapM fun v => match v with
    | .str t => parseRat t
    | _ => do let i ← v.getInt?; pure (i : ℚ)

def parseBead (j : Json) : Except String BeadSpec := do
  let ts ← parseNats (← j.getObjVal? "ts")
  let start ← (← j.getObjVal? "start").getNat?
  let stop ← (← j.getObjVal? "end").getNat?
  let steps ← match j.getObjVal? "steps" with
    | .ok (.arr a) => a.mapM (·.getNat?)
    | _ => pure #[]
  let piPoint ← match j.getObjVal? "pi_point" with
    | .ok (.arr a) => some <$> a.mapM parseMQ
    | _ => pure none
  return ⟨ts, start, stop, steps, piPoint⟩

def parseFalse (j : Json) : Except String FalseCert := do
  let isosig ← (← j.getObjVal? "isosig").getStr?
  let sqlens ← parseSqlens (← j.getObjVal? "sqlens")
  let r ← (← j.getObjVal? "reason").getStr?
  let reason : Reason ←
    if r == "invalid" then do pure (Reason.invalid (← (← j.getObjVal? "tet").getNat?))
    else if r == "edge" then do pure (Reason.edge (← (← j.getObjVal? "edge").getNat?))
    else if r == "vertex" then do
      let v ← (← j.getObjVal? "vertex").getNat?
      let nk ← (← (← j.getObjVal? "necklace").getArr?).mapM parseBead
      pure (Reason.vertex v nk)
    else throw s!"unknown reason {r}"
  return ⟨isosig, sqlens, reason⟩

def parseTree (j : Json) : Except String (Std.HashMap Nat Nat) := do
  let mut out := {}
  for (k, v) in ← objList j do
    let some x := k.toNat? | throw s!"bad bead index {k}"
    out := out.insert x (← v.getNat?)
  return out

def parseTrue (j : Json) : Except String TrueCert := do
  let isosig ← (← j.getObjVal? "isosig").getStr?
  let sqlens ← parseSqlens (← j.getObjVal? "sqlens")
  let opts := (j.getObjVal? "options").toOption.getD (Json.mkObj [])
  let piBeads := ((opts.getObjVal? "pi_beads") >>= (·.getBool?)).toOption.getD true
  let allowRepeats := ((opts.getObjVal? "allow_repeats") >>= (·.getBool?)).toOption.getD true
  let assumptions ← match j.getObjVal? "assumptions" with
    | .ok (.arr a) => a.mapM (·.getStr?)
    | _ => pure #[]
  let mut vertices := #[]
  for v in ← (← j.getObjVal? "vertices").getArr? do
    let iv ← (← v.getObjVal? "vertex").getNat?
    let nb ← (← v.getObjVal? "n_beads").getNat?
    let mut trees := {}
    for (k, t) in ← objList (← v.getObjVal? "trees") do
      let some s := k.toNat? | throw s!"bad bead index {k}"
      trees := trees.insert s (← parseTree t)
    vertices := vertices.push (iv, nb, trees)
  return ⟨isosig, sqlens, piBeads, allowRepeats, assumptions, vertices⟩

/-- Check a certificate of either kind: `(true, message)` if it proves its
claim, `(false, reason)` otherwise. Returns also which claim it makes. -/
def checkJson (j : Json) : String × Bool × String :=
  let fmt := ((j.getObjVal? "format") >>= (·.getStr?)).toOption.getD ""
  let run (r : Except String String) : Bool × String :=
    match r with
    | .ok m => (true, m)
    | .error e => (false, e)
  if fmt == FALSE_FORMAT then
    ("False", run (parseFalse j >>= checkFalse))
  else if fmt == TRUE_FORMAT then
    ("True", run (parseTrue j >>= checkTrue))
  else ("?", (false, "unknown certificate format"))

end VLink
