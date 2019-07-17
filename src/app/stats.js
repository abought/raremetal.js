/**
 * Calculate group-based tests from score statistics.
 *
 * @module stats
 * @license MIT
 */
import { cholesky } from './linalg.js';
import mvtdstpack from './mvtdstpack.js';
import numeric from 'numeric';
import * as qfc from './qfc.js';
import { GaussKronrod } from './quadrature.js';
import { pchisq, dbeta, pnorm, qchisq, dchisq } from './rstats.js';


// Functions using WASM will be defined inside a single promise- sort of a meta-module
//   Because the webassembly code is loaded asynchronously, anything using any module method will need to be
//   resolved asynchronously as well.
const MVT_WASM_HELPERS = new Promise((resolve, reject) => {
  // The emscripten "module" doesn't return a true promise, so it can't be chained in the traditional sense.
  // This syntax is a hack that allows us to wrap the wasm module with our helper functions and access those helpers.
  try {
    mvtdstpack().then(module => {
      function makeDoubleVec(size) {
        const v = new module.DoubleVec();
        v.resize(size, NaN);
        return v;
      }

      function makeIntVec(size) {
        const v = new module.IntVec();
        v.resize(size, NaN);
        return v;
      }

      function copyToDoubleVec(arr, constructor=module.DoubleVec) {
        const v = new constructor();
        for (let i = 0; i < arr.length; i++) {
          v.push_back(arr[i]);
        }
        return v;
      }

      function pmvnorm(lower, upper, mean, sigma) {
        const n = sigma.length;
        const infin = makeIntVec(n);
        const delta = makeDoubleVec(n);
        const corrF = makeDoubleVec(n * (n - 1) / 2);
        let corr = cov2cor(sigma);

        // Populate corrF
        for (let j = 0; j < n; j++) {
          for (let i = j + 1; i < n; i++) {
            let k = j + 1 + ((i - 1) * i) / 2 - 1;
            corrF.set(k, corr[i][j]);
          }
        }

        // Calculate limits
        for (let i = 0; i < n; i++) {
          delta.set(i, 0.0);

          if (lower[i] !== Infinity && lower[i] !== -Infinity) {
            lower[i] = (lower[i] - mean[i]) / Math.sqrt(sigma[i][i]);
          }

          if (upper[i] !== Infinity && upper[i] !== -Infinity) {
            upper[i] = (upper[i] - mean[i]) / Math.sqrt(sigma[i][i]);
          }

          if (lower[i] === -Infinity) {
            infin.set(i, 0);
          }
          if (upper[i] === Infinity) {
            infin.set(i, 1);
          }
          if (lower[i] === -Infinity && upper[i] === Infinity) {
            infin.set(i, -1);
          }
          if (lower[i] !== -Infinity && upper[i] !== Infinity) {
            infin.set(i, 2);
          }
          if (lower[i] === -Infinity) {
            lower[i] = 0;
          }
          if (upper[i] === Infinity) {
            upper[i] = 0;
          }
        }

        let inform = 0;
        let value = 0.0;
        let error = 0.0;
        const df = 0;
        const maxpts = 50000;
        const abseps = 0.001;
        const releps = 0.0;

        let sum = 0;
        for (let i = 0; i < n; i++) {
          sum += infin.get(i);
        }

        if (sum === -n) {
          inform = 0;
          value = 1.0;
        } else {
          ({ error, inform, value } = module.mvtdst(n, df, copyToDoubleVec(lower), copyToDoubleVec(upper), infin, corrF, delta, maxpts, abseps, releps));
        }

        if (inform === 3) {
          // Need to make correlation matrix positive definite
          let trial = 0;
          while (inform > 1 && trial < 100) {
            let eig = numeric.eig(corr, 100000);

            let lambdas = eig.lambda.x;
            for (let i = 0; i < n; i++) {
              if (lambdas[i] < 0) {
                lambdas[i] = 0.0;
              }
            }

            let D = numeric.diag(lambdas);
            let V = eig.E.x;
            corr = numeric.dot(numeric.dot(V, D), numeric.transpose(V));
            let corr_diag = Array(n);
            for (let i = 0; i < n; i++) {
              corr_diag[i] = corr[i][i];
            }
            let norm = numeric.dot(numeric.transpose([corr_diag]), [corr_diag]);

            for (let j = 0; j < n; j++) {
              for (let i = j + 1; i < n; i++) {
                let k = j + 1 + ((i - 1) * i) / 2 - 1;
                corrF.set(k, corr[i][j] / Math.sqrt(norm[i][j]));
              }
            }

            ({ error, inform, value } = module.mvtdst(n, df, copyToDoubleVec(lower), copyToDoubleVec(upper), infin, corrF, delta, maxpts, abseps, releps));
          }

          if (inform > 1) {
            value = -1.0;
          }
        }

        return {
          error: error,
          inform: inform,
          value: value
        };
      }

      const helper_module = {
        makeDoubleVec,
        makeIntVec,
        copyToDoubleVec,
        pmvnorm,
      };

      resolve(helper_module);
    });
  } catch (error) {
    reject(error);
  }
});

function emptyRowMatrix(nrows, ncols) {
  let m = new Array(nrows);
  for (let i = 0; i < nrows; i++) {
    m[i] = new Array(ncols).fill(NaN);
  }
  return m;
}

function cov2cor(sigma) {
  const corr = emptyRowMatrix(sigma.length, sigma[0].length);
  for (let i = 0; i < sigma.length; i++) {
    for (let j = i; j < sigma[0].length; j++) {
      if (i === j) {
        corr[i][j] = 1.0;
      } else {
        let v = sigma[i][j] / (Math.sqrt(sigma[i][i]) * Math.sqrt(sigma[j][j]));
        corr[i][j] = v;
        corr[j][i] = v;
      }
    }
  }
  return corr
}

function get_conditional_dist(scores, cov, comb) {
  const result = new Array(2).fill(0.0);
  const mu2 = [];
  const dim = comb.length - 1;
  const sub_cov = emptyRowMatrix(dim, dim);

  for (let i = 0; i < dim; i++) {
    let idx1 = comb[i + 1];
    mu2[i] = scores[idx1];
    for (let j = 0; j < dim; j++) {
      let idx2 = comb[j + 1];
      sub_cov[i][j] = cov[idx1][idx2];
    }
  }

  const inv = numeric.inv(sub_cov);
  const sigma12 = new Array(dim).fill(NaN);
  for (let i = 0; i < dim; i++) {
    let idx1 = comb[0];
    let idx2 = comb[i + 1];
    sigma12[i] = cov[idx1][idx2];
  }

  const tmp = new Array(dim).fill(0.0);
  for (let i = 0; i < dim; i++) {
    tmp[i] += numeric.dot(sigma12, inv[i]);
  }

  result[0] = numeric.dot(tmp, mu2);
  result[1] = 1.0 - numeric.dot(tmp, sigma12);

  if (result[1] < 0) {
    result[1] = Math.abs(result[1]);
  }

  return result;
}

/**
 * Calculates MVT p-value directly from scores/covariance and maximum test statistic.
 * TODO: ask Shaung or Goncalo where this comes from?
 * @param scores
 * @param cov_t
 * @param t_max
 * @return {*|number}
 */
function calculate_mvt_pvalue(scores, cov_t, t_max) {
  let pvalue = 0.0;
  const dim = scores.length;
  let chisq = t_max * t_max;
  let jointProbHash = {};

  if (dim === 1) {
    pvalue = pchisq(chisq, 1, 0, 0);
    return pvalue;
  }

  let uni = pchisq(chisq, 1, 0, 0);
  pvalue += dim * uni;
  let indx = [];
  let alpha = [...Array(dim).keys()]; // 0, 1, 2, 3... dim
  for (let r = 2; r <= dim; r++) {
    let j = r;
    let k = r;
    let comb = [];
    let par = [];

    for (let twk = j; twk <= k; twk++) {
      let r = twk;
      let done = true;
      for (let iwk = 0; iwk < r; iwk++) {
        indx.push(iwk);
      }

      while (done) {
        done = false;
        for (let owk = 0; owk < r; owk++) {
          comb.push(alpha[indx[owk]]);
        }

        par = get_conditional_dist(scores, cov_t, comb);
        let chisq, condProb, prob;
        if (par[1] === 0.0) {
          condProb = 0.0;
        }
        else {
          chisq = (t_max - par[0]) * (t_max - par[0]) / par[1];
          if (chisq < 0) {
            chisq = -chisq;
          }
          condProb = pchisq(chisq, 1, 0, 0);
        }

        let hashKey = "";
        if (r === 2) {
          hashKey += comb[0];
          hashKey += comb[1];
          prob = condProb * uni;
          jointProbHash[hashKey] = prob;
          hashKey = "";
        }
        else {
          for (let i = 1; i < r; i++) {
            hashKey += comb[i];
          }

          prob = jointProbHash[hashKey];
          prob *= condProb;
          let newKey = "";
          newKey += comb[0];
          newKey += hashKey;
          jointProbHash[newKey] = prob;
          hashKey = "";
        }

        pvalue -= prob;
        comb = [];
        for (let iwk = r-1; iwk >= 0; iwk--) {
          if (indx[iwk] <= (dim-1) - (r-iwk)) {
            indx[iwk]++;
            for (let swk = iwk + 1; swk < r; swk++) {
              indx[swk] = indx[swk-1] + 1;
            }
            iwk = -1;
            done = true;
          }
        }
      }
      indx = [];
    }
  }
  return pvalue;
}

/**
 * Base class for all aggregation tests.
 */
class AggregationTest {
  constructor() {
    this.label = '';
    this.key = '';

    this.requiresMaf = false;
  }

  run(u, v, w, mafs) { // todo update docstrings and call sigs
    throw new Error("Method must be implemented in a subclass");
  }
}

/**
 * Standard burden test that collapses rare variants into a total count of rare alleles observed per sample
 * in a group (e.g. gene). <p>
 *
 * See {@link https://genome.sph.umich.edu/wiki/RAREMETAL_METHOD#BURDEN_META_ANALYSIS|our wiki page} for more information.
 * Also see the {@link https://www.ncbi.nlm.nih.gov/pubmed/19810025|paper} describing the method.
 *
 * @extends AggregationTest
 */
class ZegginiBurdenTest extends AggregationTest {
  constructor() {
    super(...arguments);
    this.key = 'burden';
    this.label = 'Burden';
  }

  /**
   * Default weight function for burden test. All variants weighted equally. Only requires the number of variants
   * since they are all given the same weight value.
   * @param n {number} Number of variants.
   * @return {number[]} An array of weights, one per variant.
   */
  static weights(n) {
    return new Array(n).fill(1 / n);
  }

  /**
   * Calculate burden test from vector of score statistics and variances.
   *
   * @param {Number[]} u Vector of score statistics (length m, number of variants)
   * @param {Number[]} v Covariance matrix of score statistics
   * @param {Number[]} w Weight vector (length m, number of variants)
   * @return {Number[]} Burden test statistic z and p-value
   */
  run(u, v, w) {
    for (let e of [u, v]) {
      if (!Array.isArray(e) || !e.length) {
        throw 'Please provide all required arrays';
      }
    }

    if (!(u.length === v.length)) {
      throw 'u and v must be same length';
    }

    if (w != null) {
      if (w.length !== u.length) {
        throw 'w vector must be same length as score vector u';
      }
    }
    else {
      w = ZegginiBurdenTest.weights(u.length);
    }

    // This is taken from:
    // https://genome.sph.umich.edu/wiki/RAREMETAL_METHOD#BURDEN_META_ANALYSIS
    let over = numeric.dot(w, u);
    let under = Math.sqrt(numeric.dot(numeric.dot(w, v), w));
    let z = over / under;

    // The -Math.abs(z) is because pnorm returns the lower tail probability from the normal dist
    // The * 2 is for a two-sided p-value.
    let p = pnorm(-Math.abs(z), 0, 1) * 2;
    return [z, p];
  }
}

function _vt(maf_cutoffs, u, v, mafs) {
  // Calculate score statistic and cov weight matrix for each MAF cutoff.
  const cov_weight = emptyRowMatrix(maf_cutoffs.length, u.length);
  let t_max = -Infinity;
  const scores = Array(maf_cutoffs.length).fill(0.0);
  maf_cutoffs.map((m, i) => {
    // Weight is 1 if MAF < cutoff, otherwise 0.
    let w = mafs.map(maf => maf <= m ? 1 : 0);
    cov_weight[i] = w;

    // Calculate burden t-statistic for this maf cutoff
    let numer = numeric.dot(w, u);
    let denom = numeric.dot(numeric.dot(w, v), w);
    let t_stat = Math.abs(numer / Math.sqrt(denom));
    scores[i] = t_stat;
    if (t_stat > t_max) { t_max = t_stat; }
  });

  // Did we calculate any valid scores?
  if (Math.max(...scores) === 0.0) { throw 'No scores were able to be calculated for this group'; }

  // Calculate covariance matrix
  const cov_u = numeric.dot(numeric.dot(cov_weight, v), numeric.transpose(cov_weight));
  const cov_t = cov2cor(cov_u);

  return [scores, cov_t, t_max];
}

/**
 * Variable threshold test (VT). <p>
 */
class VTTest extends AggregationTest {
  constructor() {
    super(...arguments);
    this.label = 'Variable Threshold';
    this.key = 'vt';
    this.requiresMaf = true;
    this._method = 'auto';
  }

  /**
   * This code corresponds roughly to: https://github.com/statgen/raremetal/blob/2c82cfc5710dbd9fd56ef67a7ca5f74772d4e70d/raremetal/src/Meta.cpp#L3456
   * @param u
   * @param v
   * @param w This parameter is ignored for VT. Weights are calculated automatically from mafs.
   * @param mafs
   * @return Promise
   */
  run(u, v, w, mafs) {
    // Uses wasm, returns a promise
    if (w != null) {
      throw 'w vector is not accepted in with VT test';
    }

    // Figure out MAF cutoffs. This tries every possible MAF cutoff given a list of all MAFs.
    let maf_cutoffs = [];
    const sorted_mafs = [...mafs].sort();
    for (let i = 0; i < mafs.length; i++) {
      if (sorted_mafs[i] > maf_cutoffs.slice(-1)) {
        maf_cutoffs.push(sorted_mafs[i]);
      }
    }

    // Try calculating scores/t-stat covariance the first time (may need refinement later).
    let [scores, cov_t, t_max] = _vt(maf_cutoffs, u, v, mafs);
    const lower = new Array(maf_cutoffs.length).fill(-t_max);
    const upper = new Array(maf_cutoffs.length).fill(t_max);
    const mean = new Array(maf_cutoffs.length).fill(0);

    return MVT_WASM_HELPERS.then(module => {
      let result = module.pmvnorm(lower, upper, mean, cov_t);

      let pvalue;
      if (result.value === -1.0) {
        throw 'Error: correlation matrix is not positive semi-definite';
      } else if (result.value === 1.0) {
        // Use Shuang's algorithm
        if (maf_cutoffs.length > 20) {
          maf_cutoffs = maf_cutoffs.slice(-20);
          let [scores, cov_t, t_max] = _vt(maf_cutoffs, u, v, mafs);
          pvalue = calculate_mvt_pvalue(scores, cov_t, t_max);
        } else {
          pvalue = calculate_mvt_pvalue(scores, cov_t, t_max);
        }
      } else {
        pvalue = 1.0 - result.value;
      }

      if (pvalue > 1.0) {
        pvalue = 1.0;
      }

      return [t_max, pvalue];
    });
  }
}

/**
 * Sequence kernel association test (SKAT). <p>
 *
 * See the {@link https://www.cell.com/ajhg/fulltext/S0002-9297%2811%2900222-9|original paper} for details on the
 * method, and {@link https://genome.sph.umich.edu/wiki/RAREMETAL_METHOD#SKAT_META_ANALYSIS|our wiki} for information
 * on how the test is calculated using scores/covariances. <p>
 *
 * @extends AggregationTest
 */
class SkatTest extends AggregationTest {
  constructor() {
    super(...arguments);
    this.label = 'SKAT';
    this.key = 'skat';
    this.requiresMaf = true;

    /**
     * Skat test method. Only used for dev/testing.
     * Should not be set by user.
     * @private
     * @type {string}
     */
    this._method = 'auto';
  }

  /**
   * Calculate typical SKAT weights using beta density function.
   *
   * @function
   * @param mafs {number[]} Array of minor allele frequencies.
   * @param a {number} alpha defaults to 1.
   * @param b {number} beta defaults to 25.
   */
  static weights(mafs, a = 1, b = 25) {
    let weights = Array(mafs.length).fill(null);
    for (let i = 0; i < mafs.length; i++) {
      let w = dbeta(mafs[i], a, b, false);
      w *= w;
      weights[i] = w;
    }
    return weights;
  }

  /**
   * Calculate SKAT test. <p>
   *
   * The distribution function of the SKAT test statistic is evaluated using Davies' method by default.
   * In the special case where there is only 1 lambda, the Liu moment matching approximation method is used. <p>
   *
   * @function
   * @param {Number[]} u Vector of score statistics (length m, number of variants).
   * @param {Number[]} v Covariance matrix of score statistics (m x m).
   * @param {Number[]} w Weight vector (length m, number of variants). If weights are not provided, they will
   *  be calculated using the default weights() method of this object.
   * @param {Number[]} mafs A vector of minor allele frequencies. These will be used to calculate weights if
   *  they were not provided.
   * @return {Number[]} SKAT p-value.
   */
  run(u, v, w, mafs) {
    // Calculate weights (if necessary)
    if (w === undefined || w === null) {
      w = SkatTest.weights(mafs);
    }

    // Calculate Q
    let q = numeric.dot(numeric.dot(u,numeric.diag(w)),u);

    /**
     * Code to calculate eigenvalues from V^(1/2) * W * V^(1/2)
     * This first decomposes V = U * S * U' (SVD on symmetric normal matrix results in this, instead of U * S * V').
     * If we take sqrt(S), then U * sqrt(S) * U' is a square root of the original matrix V. For a diagonal matrix,
     * sqrt(S) is just the sqrt(s_i) of each individual diagonal element.
     * Then we just take the dot product of (U * sqrt(S) * U') * W * (U * sqrt(S) * U'), which is V^(1/2) * W * V^(1/2).
     * Finally we compute SVD of that matrix, and take the singular values as the eigenvalues.
     */
    let lambdas;
    try {
      let svd = numeric.svd(v);
      let sqrtS = numeric.sqrt(svd.S);
      let uT = numeric.transpose(svd.U);
      let eigenRhs = numeric.dot(numeric.dot(svd.U, numeric.diag(sqrtS)), uT);
      let eigenLhs = numeric.dot(eigenRhs, numeric.diag(w));
      let eigen = numeric.dot(eigenLhs, eigenRhs);
      let finalSvd = numeric.svd(eigen);
      lambdas = numeric.abs(finalSvd.S);
    } catch(error) {
      console.log(error);
      return [NaN, NaN];
    }

    if (numeric.sum(lambdas) < 0.0000000001) {
      console.error("Sum of lambda values for SKAT test is essentially zero");
      return [NaN, NaN];
    }

    // P-value method
    if (this._method === 'liu') {
      // Only for debug purposes
      return _skatLiu(lambdas, q);
    }
    else if (this._method === 'davies') {
      return _skatDavies(lambdas, q);
    }
    else if (this._method === 'auto') {
      if (lambdas.length === 1) {
        // Davies method does not support 1 lambda
        // This is what raremetal does
        return _skatLiu(lambdas, q);
      }
      else {
        let daviesResult = _skatDavies(lambdas, q);
        if (isNaN(daviesResult[1])) {
          // Davies' method could not converge. Use R-SKAT's approach instead.
          return _skatLiu(lambdas, q);
        } else {
          return daviesResult;
        }
      }
    }
    else {
      throw new Error(`Skat method ${this._method} not implemented`);
    }
  }
}

/**
 * Calculate SKAT p-value using Davies method.
 * @function
 * @param lambdas Eigenvalues of sqrtV * W * sqrtV.
 * @param qstat SKAT test statistic U.T * W * U.
 * @return {Number[]} Array of [Q statistic, p-value].
 * @private
 */
function _skatDavies(lambdas, qstat) {
  /**
   * lambdas - coefficient of jth chi-squared variable
   * nc1 - non-centrality parameters
   * n1 - degrees of freedom
   * n - number of chi-squared variables
   * sigma - coefficient of standard normal variable
   * qstat - point at which cdf is to be evaluated (this is SKAT Q stat usually)
   * lim1 - max number of terms in integration
   * acc - maximum error
   * trace - array into which the following is stored:
   *   trace[0]	absolute sum
   *   trace[1]	total number of integration terms
   *   trace[2]	number of integrations
   *   trace[3]	integration interval in final integration
   *   trace[4]	truncation point in initial integration
   *   trace[5]	s.d. of initial convergence factor
   *   trace[6]	cycles to locate integration parameters
   * ifault - array into which the following fault codes are stored:
   *   0 normal operation
   *   1 required accuracy not achieved
   *   2 round-off error possibly significant
   *   3 invalid parameters
   *   4 unable to locate integration parameters
   *   5 out of memory
   * res - store final value into this variable
   */
  let n = lambdas.length;
  let nc1 = Array(n).fill(0);
  let n1 = Array(n).fill(1);
  let sigma = 0.0;
  let lim1 = 10000;
  let acc = 0.0001;
  let res = qfc.qf(lambdas, nc1, n1, n, sigma, qstat, lim1, acc);
  let qfval = res[0];
  let pval = 1.0 - qfval;

  let converged = (res[1] === 0) && (pval > 0) && (pval <= 1);
  if (!converged) {
    pval = NaN;
  }

  return [qstat, pval];
}

/**
 * Calculate SKAT p-value using Liu method.
 * @param lambdas Eigenvalues of sqrtV * W * sqrtV.
 * @param qstat SKAT test statistic U.T * W * U.
 * @return {Number[]} [qstat, pvalue]
 * @private
 */
function _skatLiu(lambdas, qstat) {
  let n = lambdas.length;
  let [c1, c2, c3, c4] = Array(4).fill(0.0);
  for (let i = 0; i < n; i++) {
    let ilambda = lambdas[i];
    c1 += ilambda;
    c2 += ilambda * ilambda;
    c3 += ilambda * ilambda * ilambda;
    c4 += ilambda * ilambda * ilambda * ilambda;
  }

  let s1 = c3 / Math.sqrt(c2 * c2 * c2);
  let s2 = c4 / (c2 * c2);
  let muQ = c1;
  let sigmaQ = Math.sqrt(2.0 * c2);
  let tStar = (qstat - muQ) / sigmaQ;

  let delta, l, a;
  if (s1 * s1 > s2) {
    a = 1.0 / (s1 - Math.sqrt(s1 * s1 - s2));
    delta = s1 * a * a * a - a * a;
    l = a * a - 2.0 * delta;
  } else {
    a = 1.0 / s1;
    delta = 0.0;
    l = c2 * c2 * c2 / (c3 * c3);
  }

  let muX = l + delta;
  let sigmaX = Math.sqrt(2.0) * a;
  let qNew = tStar * sigmaX + muX;
  let p;

  if (delta === 0) {
    p = pchisq(qNew,l,0,0);
  } else {
    // Non-central chi-squared
    p = pchisq(qNew,l,delta,0,0);
  }

  return [qstat, p];
}

function getEigen(m) {
  const lambdas = numeric.eig(m, 1000000).lambda.x.sort((a, b) => a - b);
  const n = lambdas.length;
  let numNonZero = 0;
  let sumNonZero = 0.0;
  for (let i = 0; i < n; i++) {
    if (lambdas[i] > 0) {
      numNonZero++;
      sumNonZero += lambdas[i];
    }
  }

  if (numNonZero === 0) {
    throw new Error("All eigenvalues were 0 when calculating SKAT-O test");
  }

  const t = sumNonZero / numNonZero / 100000;
  let numKeep = n;
  for (let i = 0; i < n; i++) {
    if (lambdas[i] < t) {
      numKeep--;
    }
    else {
      break;
    }
  }

  const keep = new Array(numKeep).fill(null);
  for (let i = 0; i < numKeep; i++) {
    keep[i] = lambdas[n - 1 - i];
  }

  return keep;
}

function getMoment(lambdas) {
  let c = new Array(4).fill(NaN);
  c[0] = numeric.sum(lambdas);
  c[1] = numeric.sum(numeric.pow(lambdas, 2));
  c[2] = numeric.sum(numeric.pow(lambdas, 3));
  c[3] = numeric.sum(numeric.pow(lambdas, 4));

  const muQ = c[0];
  const sigmaQ = Math.sqrt(2 * c[1]);
  const s1 = c[2] / c[1] / Math.sqrt(c[1]);
  const s2 = c[3] / (c[1] * c[1]);

  let a, d, l;
  if (s1 * s1 > s2) {
    a = 1 / (s1 - Math.sqrt(s1 * s1 - s2));
    d = (s1 * a - a * a);
    l = a * a - 2 * d;
  }
  else {
    l = 1.0 / s2;
  }

  const varQ = sigmaQ * sigmaQ;
  const df = l;
  return {
    muQ: muQ,
    varQ: varQ,
    df: df
  }
}

function getPvalByMoment(q, m) {
  const qNorm = (q - m.muQ) / Math.sqrt(m.varQ) * Math.sqrt(2.0 * m.df) + m.df;
  return pchisq(qNorm, m.df, 0, false, false);
}

function getQvalByMoment(min_pval, m) {
  const q_org = qchisq(min_pval, m.df, 0, false, false);
  return (q_org - m.df) / Math.sqrt(2.0 * m.df) * Math.sqrt(m.varQ) + m.muQ;
}

class SkatIntegrator {
  constructor(rhos, lambda, Qs_minP, taus, MuQ, VarQ, VarZeta, Df) {
    this.rhos = rhos;
    this.lambda = lambda;
    this.Qs_minP = Qs_minP;
    this.taus = taus;
    this.MuQ = MuQ;
    this.VarQ = VarQ;
    this.VarZeta = VarZeta;
    this.Df = Df;
  }

  static pvalueDavies(q, lambdas) {
    let n = lambdas.length;
    let nc1 = Array(n).fill(0);
    let n1 = Array(n).fill(1);
    let sigma = 0.0;
    let lim1 = 10000;
    let acc = 0.0001;
    let res = qfc.qf(lambdas, nc1, n1, n, sigma, q, lim1, acc);
    let qfval = res[0];
    let fault = res[1];
    let pvalue = 1.0 - qfval;

    if (pvalue > 1.0) {
      pvalue = 1.0;
    }

    if (fault) {
      pvalue = -1.0;
    }

    return pvalue;
  }

  static pvalueLiu(q, lambdas) {
    let n = lambdas.length;
    let [c1, c2, c3, c4] = Array(4).fill(0.0);
    for (let i = 0; i < n; i++) {
      let ilambda = lambdas[i];
      c1 += ilambda;
      c2 += ilambda * ilambda;
      c3 += ilambda * ilambda * ilambda;
      c4 += ilambda * ilambda * ilambda * ilambda;
    }

    let s1 = c3 / Math.sqrt(c2 * c2 * c2);
    let s2 = c4 / (c2 * c2);
    let muQ = c1;
    let sigmaQ = Math.sqrt(2.0 * c2);
    let tStar = (q - muQ) / sigmaQ;

    let delta, l, a;
    if (s1 * s1 > s2) {
      a = 1.0 / (s1 - Math.sqrt(s1 * s1 - s2));
      delta = s1 * a * a * a - a * a;
      l = a * a - 2.0 * delta;
    } else {
      a = 1.0 / s1;
      delta = 0.0;
      l = c2 * c2 * c2 / (c3 * c3);
    }

    let muX = l + delta;
    let sigmaX = Math.sqrt(2.0) * a;
    let qNew = tStar * sigmaX + muX;

    if (qNew < 0) { return 1; }

    let p;
    if (delta === 0) {
      p = pchisq(qNew,l,0,0);
    } else {
      // Non-central chi-squared
      p = pchisq(qNew,l,delta,0,0);
    }

    return p;
  }

  integrandDavies(x) {
    let kappa = Number.MAX_VALUE;
    const nRho = this.rhos.length;
    for (let i = 0; i < nRho; i++) {
      let v = (this.Qs_minP[i] - this.taus[i] * x) / (1.0 - this.rhos[i]);
      if (i === 0) {
        kappa = v;
      }
      if (v < kappa) {
        kappa = v;
      }
    }
    let temp;
    if (kappa > numeric.sum(this.lambda) * 10000) {
      temp = 0.0;
    }
    else {
      let Q = (kappa - this.MuQ) * Math.sqrt(this.VarQ - this.VarZeta) / Math.sqrt(this.VarQ) + this.MuQ;
      temp = SkatIntegrator.pvalueDavies(Q, this.lambda);
      if (temp <= 0.0 || temp === 1.0) {
        temp = SkatIntegrator.pvalueLiu(Q, this.lambda);
      }
    }
    let final = (1.0 - temp) * dchisq(x, 1);
    //console.log("integrandDavies: ", x, temp, final);
    return final;
  }

  integrandLiu(x) {
    let kappa = Number.MAX_VALUE;
    const nRho = this.rhos.length;
    for (let i = 0; i < nRho; i++) {
      let v = (this.Qs_minP[i] - this.taus[i] * x) / (1.0 - this.rhos[i]);
      if (v < kappa) {
        kappa = v;
      }
    }
    let Q = (kappa - this.MuQ) / Math.sqrt(this.VarQ) * Math.sqrt(2.0 * this.Df) + this.Df;

    let ret;
    if (Q <= 0) {
      ret = 0;
    }
    else {
      ret = pchisq(Q, this.Df) * dchisq(x, 1);
    }

    return ret;
  }

  skatOptimalIntegral() {
    const integ = new GaussKronrod(21, 15);

    // Try integrating Davies first
    let result;
    try {
      result = integ.integrate(this.integrandDavies.bind(this), 0, 40);
    }
    catch (e1) {
      try {
        result = integ.integrate(this.integrandLiu.bind(this), 0, 40);
      }
      catch (e2) {
        console.error("Could not integrate Davies or Liu integrands (SKAT-O)");
        throw e2;
      }
    }

    return result[0];
  }
}

/**
 * Optimal sequence kernel association test (SKAT). <p>
 *
 * The following papers detail the method:
 *
 * Original SKAT optimal test paper, utilizing genotypes instead of covariance matrices: https://doi.org/10.1016/j.ajhg.2012.06.007
 * Meta-analysis of SKAT optimal test, and use of covariance matrices: https://doi.org/10.1016/j.ajhg.2013.05.010
 *
 * @extends AggregationTest
 */
class SkatOptimalTest extends AggregationTest {
  constructor() {
    super(...arguments);
    this.label = 'SKAT Optimal';
    this.key = 'skat-o';
    this.requiresMaf = true;

    /**
     * Skat test method. Only used for dev/testing.
     * Should not be set by user.
     * @private
     * @type {string}
     */
    this._method = 'auto';
  }

  /**
   * Calculate typical SKAT weights using beta density function.
   *
   * @function
   * @param mafs {number[]} Array of minor allele frequencies.
   * @param a {number} alpha defaults to 1.
   * @param b {number} beta defaults to 25.
   */
  static weights(mafs, a = 1, b = 25) {
    let weights = Array(mafs.length).fill(null);
    for (let i = 0; i < mafs.length; i++) {
      let w = dbeta(mafs[i], a, b, false);
      //w *= w;
      weights[i] = w;
    }
    return weights;
  }

  /**
   * Calculate optimal SKAT test. <p>
   *
   * @function
   * @param {Number[]} u Vector of score statistics (length m, number of variants).
   * @param {Number[]} v Covariance matrix of score statistics (m x m).
   * @param {Number[]} w Weight vector (length m, number of variants). If weights are not provided, they will
   *  be calculated using the default weights() method of this object.
   * @param {Number[]} mafs A vector of minor allele frequencies. These will be used to calculate weights if
   *  they were not provided.
   * @param {Number[]} rhos A vector of rho values, representing the weighting between burden and SKAT statistics.
   * @return {Number[]} SKAT p-value.
   */
  run(u, v, w, mafs, rhos) {
    const { dot, sum, mul, div, sub, rep, pow, diag } = numeric;
    const t = numeric.transpose;

    if (u.length === 1) {
      // rvtest
      return new SkatTest().run(u, v, w, mafs);
    }

    // Calculate weights (if necessary)
    if (w === undefined || w === null) {
      w = SkatOptimalTest.weights(mafs);
    }

    const nVar = u.length; // number of variants
    w = diag(w); // diagonal matrix
    u = t([u]); // column vector

    // Setup rho values
    if (!rhos) {
      rhos = [];
      for (let i = 0; i <= 10; i++) {
        let v = i / 10;
        if (v > 0.999) {
          // rvtests does this to avoid rank deficiency
          v = 0.999;
        }
        rhos.push(v);
      }
    }
    const nRhos = rhos.length;
    // MetaSKAT optimal.mod rho values
    //const rhos = [0, 0.01, 0.04, 0.09, 0.25, 0.5, 0.999];
    //const nRhos = rhos.length;

    // Calculate rho matrices (1-rho)*I + rho*1*1'
    // [ 1   rho rho ]
    // [ rho 1   rho ]
    // [ rho rho 1   ]
    const Rp = new Array(nRhos).fill(null);
    for (let i = 0; i < nRhos; i++) {
      let r = rep([nVar, nVar], rhos[i]);
      for (let j = 0; j < r.length; j++) {
        r[j][j] = 1.0;
      }
      Rp[i] = r;
    }

    // Calculate Q statistics, where Q = U' * W * R(rho) * W * U
    // U is the score statistic vector, W is the diagonal weight matrix for each variant
    // R(rho) is a matrix for each rho value that reflects weighting between burden & SKAT
    const Qs = [];
    for (let i = 0; i < nRhos; i++) {
      Qs[i] = dot(t(u), dot(w, dot(Rp[i], dot(w, u))))[0][0];
      Qs[i] = Qs[i] / 2.0; // SKAT R package divides by 2
    }

    // Calculate lambdas (eigenvalues of W * IOTA * W.) In the paper, IOTA is the covariance matrix divided by
    // the phenotypic variance sigma^2.
    const lambdas = new Array(nRhos).fill(null);
    const phi = div(dot(w, dot(v, w)), 2); // https://git.io/fjwqF
    for (let i = 0; i < nRhos; i++) {
      let L = cholesky(Rp[i]);
      let phi_rho = dot(t(L), dot(phi, L));
      try {
        lambdas[i] = getEigen(phi_rho);
      }
      catch (error) {
        console.error(error.message);
        return [NaN, NaN];
      }
    }

    // Calculate moments
    const moments = new Array(nRhos).fill(null);
    for (let i = 0; i < nRhos; i++) {
      moments[i] = getMoment(lambdas[i]);
    }

    // Calculate p-values for each rho
    const pvals = new Array(nRhos).fill(null);
    for (let i = 0; i < nRhos; i++) {
      pvals[i] = getPvalByMoment(Qs[i], moments[i]);
    }

    // Calculate minimum p-value across all rho values
    let minP = pvals[0];
    let minIndex = 0;
    for (let i = 1; i < nRhos; i++) {
      if (pvals[i] < minP) {
        minP = pvals[i];
        minIndex = i;
      }
    }
    //const rho = rhos[minIndex];
    const Q = Qs[minIndex];

    // Calculate minimum Q(p)
    const Qs_minP = new Array(nRhos).fill(null);
    for (let i = 0; i < nRhos; i++) {
      Qs_minP[i] = getQvalByMoment(minP, moments[i]);
    }

    // Calculate parameters needed for Z'(I-M)Z part
    const Z11 = dot(phi, rep([nVar, 1], 1));
    const ZZ = phi;
    const ZMZ = div(dot(Z11, t(Z11)),sum(ZZ));
    const ZIMZ = sub(ZZ,ZMZ);
    let lambda;
    try {
      lambda = getEigen(ZIMZ);
    }
    catch (error) {
      console.error(error.message);
      return [NaN, NaN];
    }
    const varZeta = 4 * sum(mul(ZMZ, ZIMZ));
    const muQ = sum(lambda);
    const varQ = 2.0 * sum(pow(lambda, 2)) + varZeta;
    const kerQ = 12.0 * sum(pow(lambda, 4)) / ((sum(pow(lambda, 2))) ** 2);
    const dF = 12.0 / kerQ;

    // Calculate taus
    const z_mean = sum(ZZ) / (nVar ** 2);
    const tau1 = sum(dot(ZZ, ZZ)) / (nVar ** 2) / z_mean;
    const taus = new Array(nRhos).fill(null);
    for (let i = 0; i < nRhos; i++) {
      taus[i] = (nVar * nVar) * rhos[i] * z_mean + tau1 * (1 - rhos[i]);
    }

    // Calculate final p-value
    if (new Set([rhos.length, Qs_minP.length, taus.length]).size > 1) {
      throw "Parameter arrays for SKAT integration must all be the same length";
    }

    const integrator = new SkatIntegrator(
      rhos,
      lambda,
      Qs_minP,
      taus,
      muQ,
      varQ,
      varZeta,
      dF
    );
    let pvalue = 1 - integrator.skatOptimalIntegral();

    // Check SKAT p-value
    const multi = (nRhos < 3) ? 2 : 3;
    if (nRhos) {
      if (pvalue <= 0) {
        let p = minP * multi;
        if (pvalue < p) {
          pvalue = p;
        }
      }
    }
    if (pvalue === 0.0) {
      pvalue = pvals[0];
      for (let i = 1; i < nRhos; i++) {
        if (pvals[i] > 0 && pvals[i] < pvalue) {
          pvalue = pvals[i];
        }
      }
    }
    return [Q, pvalue];
  }
}

export { // for unit testing only
  AggregationTest as _AggregationTest,
  get_conditional_dist as _get_conditional_dist
};
export { SkatTest, SkatOptimalTest, ZegginiBurdenTest, VTTest,
  MVT_WASM_HELPERS, calculate_mvt_pvalue, _skatDavies, _skatLiu,

};
