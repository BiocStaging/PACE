// fp_no_contract.hpp -- forbid fused multiply-add (FMA) contraction in the code
// that follows, so distances are computed as separate rounded operations, the
// same way R and dbscan::frNN compute them on this platform.
//
// Include this FIRST in any translation unit that computes distances or kernel
// weights that must match R bit for bit. The distance code also keeps each
// multiply and add in its own statement, which already prevents contraction
// under clang's default (-ffp-contract=on); the pragmas cover compilers that
// contract across statements (gcc's default for GNU dialects on FMA targets).
#ifndef PACE_FP_NO_CONTRACT_HPP
#define PACE_FP_NO_CONTRACT_HPP

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#endif  // PACE_FP_NO_CONTRACT_HPP
