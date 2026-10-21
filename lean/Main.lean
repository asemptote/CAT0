import VLink

/-!
`vlinkcheck FILE.jsonl`: check every certificate in FILE (one JSON object per
line, or a `certificate` field of one). Prints one line per certificate:
`ACCEPT <claim> <message>` or `REJECT <claim> <reason>`.
-/

open Lean (Json)

def main (args : List String) : IO UInt32 := do
  let [path] := args | do
    IO.eprintln "usage: vlinkcheck FILE.jsonl"
    return 2
  let mut nAcc := 0
  let mut nRej := 0
  for line in (← IO.FS.lines path) do
    if line.trimAscii.toString.isEmpty then continue
    match Json.parse line with
    | .error e =>
      IO.println s!"REJECT ? malformed JSON: {e}"
      nRej := nRej + 1
    | .ok j =>
      let j := (j.getObjVal? "certificate").toOption.getD j
      let (claim, ok, msg) := VLink.checkJson j
      IO.println s!"{if ok then "ACCEPT" else "REJECT"} {claim} {msg}"
      if ok then nAcc := nAcc + 1 else nRej := nRej + 1
  IO.eprintln s!"{nAcc} accepted, {nRej} rejected"
  return 0
