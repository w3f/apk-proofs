## Recursive validator set validation

See comments in recursive.rs

> cargo run --release --features "parallel print-trace" --example recursive 8 2

will produce a fancy output like

```
Running a chain with 2^8-1 validators for 2 eras. To change the values run with '--example recursive LOG_N N_ERAS'

Setup: max validator set size = 2^8-1

Start:   Generating URS to support 2^8-1 signers
··Start:   Computing 766 scalars powers
··End:     Computing 766 scalars powers ............................................67.700µs
··Start:   766-scalar mul in G1
··End:     766-scalar mul in G1 ....................................................62.327ms
··Start:   2-scalar mul in G1
··End:     2-scalar mul in G1 ......................................................26.416ms
End:     Generating URS to support 2^8-1 signers ...................................92.094ms

Genesis: validator set size = 255, quorum = 171

Start:   Computing commitment to the set of initial 255 validators
End:     Computing commitment to the set of initial 255 validators .................92.884ms

Era 1

Start:   Each (honest) validators computes the commitment to the new validator set of size 255 and signs the commitment
End:     Each (honest) validators computes the commitment to the new validator set of size 255 and signs the commitment 637.104ms

Start:   Helper aggregated 234 individual signatures on the same commitment and generates accountable light client proof of them
End:     Helper aggregated 234 individual signatures on the same commitment and generates accountable light client proof of them 241.735ms

Start:   Light client verifies light client proof for 234 signers
··Start:   apk proof verification
····Start:   linear accountability check
····End:     linear accountability check ...........................................86.900µs
····Start:   KZG check
······Start:   linearization polynomial commitment
······End:     linearization polynomial commitment .................................1.849ms
······Start:   aggregate evaluation claims in zeta
······End:     aggregate evaluation claims in zeta .................................837.200µs
······Start:   batched KZG openning
······End:     batched KZG openning ................................................8.882ms
······Start:   lazy subgroup check
······End:     lazy subgroup check .................................................595.100µs
····End:     KZG check .............................................................12.538ms
··End:     apk proof verification ..................................................14.160ms
··Start:   aggregate BLS signature verification
··End:     aggregate BLS signature verification ....................................3.154ms
End:     Light client verifies light client proof for 234 signers ..................17.491ms

Era 2

Start:   Each (honest) validators computes the commitment to the new validator set of size 255 and signs the commitment
End:     Each (honest) validators computes the commitment to the new validator set of size 255 and signs the commitment 617.583ms

Start:   Helper aggregated 221 individual signatures on the same commitment and generates accountable light client proof of them
End:     Helper aggregated 221 individual signatures on the same commitment and generates accountable light client proof of them 229.157ms

Start:   Light client verifies light client proof for 221 signers
··Start:   apk proof verification
····Start:   linear accountability check
····End:     linear accountability check ...........................................87.200µs
····Start:   KZG check
······Start:   linearization polynomial commitment
······End:     linearization polynomial commitment .................................1.888ms
······Start:   aggregate evaluation claims in zeta
······End:     aggregate evaluation claims in zeta .................................902.900µs
······Start:   batched KZG openning
······End:     batched KZG openning ................................................8.688ms
······Start:   lazy subgroup check
······End:     lazy subgroup check .................................................629.400µs
····End:     KZG check .............................................................12.608ms
··End:     apk proof verification ..................................................14.196ms
··Start:   aggregate BLS signature verification
··End:     aggregate BLS signature verification ....................................3.178ms
End:     Light client verifies light client proof for 221 signers ..................17.646ms
```