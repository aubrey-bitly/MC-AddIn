/* =============================================================================
   Monte Carlo Simulator v1 — Custom Functions (functions.js)
   Plain JavaScript — no build tools required.
   Registered via CustomFunctions.associate() API.
   Namespace: MC  =>  =MC.PERT(...), =MC.SAMPLE(...), etc.
   ============================================================================= */

"use strict";

// =============================================================================
//  MATH LIBRARY (self-contained — all algorithms inline)
// =============================================================================

// ---- Log-Gamma (Lanczos approximation) ----
function logGamma(x) {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
             -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) ser += c[j] / ++y;
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

// ---- Beta continued fraction (Lentz's method) ----
function betaCF(a, b, x) {
  const MAXIT = 200, EPS = 3e-12, FPMIN = 1e-30;
  const qab = a + b, qap = a + 1, qam = a - 1;
  let c = 1, d = 1 - qab * x / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return h;
}

// ---- Regularized incomplete Beta function ----
function betaInc(a, b, x) {
  if (x < 0 || x > 1) throw new Error("betaInc: x out of [0,1]");
  if (x === 0) return 0;
  if (x === 1) return 1;
  const lbeta = logGamma(a) + logGamma(b) - logGamma(a + b);
  const front = Math.exp(Math.log(x) * a + Math.log(1 - x) * b - lbeta) / a;
  if (x < (a + 1) / (a + b + 2)) return front * betaCF(a, b, x);
  return 1 - Math.exp(Math.log(1 - x) * b + Math.log(x) * a - lbeta) / b * betaCF(b, a, 1 - x);
}

// ---- Beta inverse CDF via Newton-Raphson ----
function betaInv(p, a, b) {
  if (p <= 0) return 0;
  if (p >= 1) return 1;
  const lbeta = logGamma(a) + logGamma(b) - logGamma(a + b);
  let x = a / (a + b); // start at mean
  for (let i = 0; i < 100; i++) {
    const err = betaInc(a, b, x) - p;
    if (Math.abs(err) < 1e-12) break;
    const pdf = Math.exp((a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) - lbeta);
    if (pdf < 1e-300) break;
    x -= err / pdf;
    x = Math.max(1e-12, Math.min(1 - 1e-12, x));
  }
  return x;
}

// ---- PERT distribution — inverse CDF ----
function pertInv(p, min, likely, max, lambda) {
  lambda = (lambda === undefined || lambda === null) ? 4 : lambda;
  if (min === max) return min;
  if (p <= 0) return min;
  if (p >= 1) return max;
  const range = max - min;
  const a1 = 1 + lambda * (likely - min) / range;
  const a2 = 1 + lambda * (max - likely) / range;
  return min + range * betaInv(p, a1, a2);
}

// ---- Normal inverse CDF (Peter Acklam's rational approximation) ----
// Max error < 1.15e-9 over the full (0,1) range.
function normalInv(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  const a = [-3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,
              1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00];
  const b = [-5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,
              6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
             -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00];
  const d = [ 7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,
              3.754408661907416e+00];
  const pLow = 0.02425, pHigh = 1 - pLow;
  let x;
  if (p < pLow) {
    const q = Math.sqrt(-2 * Math.log(p));
    x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    const q = p - 0.5, r = q * q;
    x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5]) * q /
        (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    const q = Math.sqrt(-2 * Math.log(1 - p));
    x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
         ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
  return x;
}

// =============================================================================
//  CUSTOM FUNCTION IMPLEMENTATIONS
// =============================================================================

/**
 * =MC.PERT(min, likely, max, lambda, percentile)
 * Returns the value at the given percentile (0–100) of a PERT distribution.
 * @customfunction
 */
CustomFunctions.associate("PERT", function(min, likely, max, lambda, percentile) {
  if (typeof min !== "number" || typeof likely !== "number" ||
      typeof max !== "number" || typeof lambda !== "number" ||
      typeof percentile !== "number") {
    throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
      "All arguments must be numbers.");
  }
  if (min >= max) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "min must be less than max.");
  if (likely < min || likely > max) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "likely must be between min and max.");
  if (percentile < 0 || percentile > 100) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "percentile must be between 0 and 100.");
  if (lambda <= 0) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "lambda must be > 0.");
  return pertInv(percentile / 100, min, likely, max, lambda);
});

/**
 * =MC.PERTMEAN(min, likely, max, lambda)
 * Returns the mean (expected value) of the PERT distribution.
 * Formula: (min + lambda * likely + max) / (lambda + 2)
 * @customfunction
 */
CustomFunctions.associate("PERTMEAN", function(min, likely, max, lambda) {
  if (typeof min !== "number" || typeof likely !== "number" ||
      typeof max !== "number" || typeof lambda !== "number") {
    throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
      "All arguments must be numbers.");
  }
  if (lambda <= 0) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "lambda must be > 0.");
  return (min + lambda * likely + max) / (lambda + 2);
});

/**
 * =MC.PERTSD(min, likely, max, lambda)
 * Returns the standard deviation of the PERT distribution.
 * @customfunction
 */
CustomFunctions.associate("PERTSD", function(min, likely, max, lambda) {
  if (typeof min !== "number" || typeof likely !== "number" ||
      typeof max !== "number" || typeof lambda !== "number") {
    throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
      "All arguments must be numbers.");
  }
  if (min >= max) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "min must be less than max.");
  if (lambda <= 0) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "lambda must be > 0.");
  const mu = (min + lambda * likely + max) / (lambda + 2);
  const range = max - min;
  return Math.sqrt((mu - min) * (max - mu)) / Math.sqrt(lambda + 3);
});

/**
 * =MC.NORMALINV(probability)
 * Returns the inverse of the standard normal CDF.
 * Equivalent to Excel's NORM.S.INV(). Max error < 1.15e-9.
 * @customfunction
 */
CustomFunctions.associate("NORMALINV", function(p) {
  if (typeof p !== "number") {
    throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
      "Argument must be a number.");
  }
  if (p <= 0 || p >= 1) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "Probability must be strictly between 0 and 1.");
  return normalInv(p);
});

/**
 * =MC.BETAINV(probability, alpha, beta)
 * Returns the inverse of the Beta distribution CDF.
 * Finds x such that BetaCDF(x; alpha, beta) = probability.
 * @customfunction
 */
CustomFunctions.associate("BETAINV", function(p, a, b) {
  if (typeof p !== "number" || typeof a !== "number" || typeof b !== "number") {
    throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
      "All arguments must be numbers.");
  }
  if (p <= 0 || p >= 1) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "Probability must be strictly between 0 and 1.");
  if (a <= 0 || b <= 0) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "Alpha and beta must be > 0.");
  return betaInv(p, a, b);
});

/**
 * =MC.SAMPLE(min, likely, max, lambda)
 * Returns a single random PERT sample. VOLATILE — recalculates on every change.
 * For large simulations use the task pane instead.
 * @customfunction
 * @volatile
 */
CustomFunctions.associate("SAMPLE", function(min, likely, max, lambda) {
  if (typeof min !== "number" || typeof likely !== "number" ||
      typeof max !== "number" || typeof lambda !== "number") {
    throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
      "All arguments must be numbers.");
  }
  if (min >= max) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "min must be less than max.");
  if (likely < min || likely > max) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "likely must be between min and max.");
  if (lambda <= 0) throw new CustomFunctions.Error(CustomFunctions.ErrorCode.invalidValue,
    "lambda must be > 0.");
  return pertInv(Math.random(), min, likely, max, lambda);
});
