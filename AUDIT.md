## R function coverage audit

### Coverage summary

| Category | Ported | Total | % |
|---|---|---|---|
| Exported R functions | 7 | 13 | 53.8% |
| Internal helpers (reachable) | 0 | 0 | 0.0% |

_Python package exposes 14 unique names._

### Exported R functions

| R function | Python equivalent | Status |
|---|---|---|
| `create_differential_topology` | `create_differential_topology` | ✅ ported |
| `differentiationTest` | `differentiationTest` | ✅ ported |
| `fateSelectionTest` | `—` | ❌ MISSING |
| `fateSelectionTest_multipleSamples` | `—` | ❌ MISSING |
| `imbalance_score` | `imbalance_score` | ✅ ported |
| `merge_sds` | `merge_sds` | ✅ ported |
| `nLineages` | `—` | ❌ MISSING |
| `progressionTest` | `progressionTest` | ✅ ported |
| `progressionTest_multipleSamples` | `—` | ❌ MISSING |
| `slingshot_conditions` | `—` | ❌ MISSING |
| `topologyTest` | `topologyTest` | ✅ ported |
| `topologyTest_multipleSamples` | `—` | ❌ MISSING |
| `weights_from_pst` | `weights_from_pst` | ✅ ported |

### Internal helpers reachable from exports

| R helper | File | Python equivalent | Status |
|---|---|---|---|
