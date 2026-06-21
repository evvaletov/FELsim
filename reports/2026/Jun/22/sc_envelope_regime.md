# Prescribed-envelope SC: regime of validity + multi-cavity sig_env

Author: Eremey Valetov, 2026-06-20.

The optional prescribed-envelope SC uses a fixed per-cell sigma, so it is accurate only while the beam stays near that envelope.

## Charge scan (SLAC linac)
| charge (nC) | prescribed/self-consistent eps_n rel diff |
|---:|---:|
| 0.001 | 0.0% |
| 0.01 | 0.2% |
| 0.03 | 1.5% |
| 0.1 | 15.7% |
| 0.3 | 108.4% |
| 1.0 | 551.4% |

At low charge (0.001 nC, weak SC, stable beam) the prescribed envelope agrees with self-consistent to 0.0%; at 1.0 nC (strong unmatched SC) it diverges to 551% because a fixed sigma over-defocuses once the beam moves. Transverse focusing keeps the beam near the envelope and extends the agreement -- that focusing-lattice validation is deferred to the real injector+linac (a designed matched FODO; a quick synthetic one was not stable).

## Multi-cavity sig_env
Two contiguous cavities build with a flat 28-cell envelope indexed by a running offset; each cavity's SC uses the correct slice (verified). Limitation: a non-cavity SC element (a quad/drift slice) between cavities is not counted by the cavity-only offset, so a fully interspersed accelerating-FODO would need the envelope to span all SC cells -- future work. update_on_track stays the default.