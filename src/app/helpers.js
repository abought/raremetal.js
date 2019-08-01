/**
 * Helper methods for running aggregation tests
 *
 * This wraps internal functionality and provides utilities for reading and writing expected API formats
 */
import numeric from "numeric";
import { REGEX_EPACTS } from "./constants";
import { _AggregationTest, SkatTest, ZegginiBurdenTest, VTTest, SkatOptimalTest } from "./stats";

const _all_tests = [ZegginiBurdenTest, SkatTest, VTTest, SkatOptimalTest];

/**
 * Look up aggregation tests by unique name.
 *
 * This is a helper for external libraries; it provides an immutable registry of all available tests.
 * TODO would be nice to get rid of this?
 *
 *
 * {key: {label: String, constructor: Object }
 * @type {{String: {label: String, constructor: function}}}
 */
const AGGREGATION_TESTS = Object.freeze(_all_tests.reduce(function (acc, constructor) {
  const inst = new constructor();  // Hack- need instance to access attributes
  acc[inst.key] = { label: inst.label, constructor: constructor };
  return acc;
}, {}));


/**
 * Helper object for reading and interpreting variant data
 */
class PortalVariantsHelper {
  constructor(variants_array) {
    this._variants = variants_array;
    this._variant_lookup = this.parsePortalVariantData(variants_array);
  }

  get data() {  // Raw unparsed data
    return this._variants;
  }

  parsePortalVariantData(variants) {
    // Read an array of variants. Parse names into position/ref/alt, and assign altFreq to MAF.
    // Return a hash keyed on variant ID for quick lookups.
    let lookup = {};
    variants.forEach(data => {
      let { variant, altFreq, pvalue, score } = data;
      let [_, chrom, pos, ref, alt, __] = variant.match(REGEX_EPACTS);  // eslint-disable-line no-unused-vars

      let effectFreq = altFreq;
      let effect = alt;

      /**
       * The variant's score statistic in the API is coded toward the alternate allele.
       * However, we want the effect coded towards the minor allele, since most rare variant tests assume
       * you are counting the rare/minor allele.
       */
      if (altFreq > 0.5) {
        /**
         * The effect allele is initially the alt allele. Since we're flipping it,
         * the "other" allele is the reference allele.
         */
        score = -score;
        effect = ref;

        // This is also now the minor allele frequency.
        effectFreq = 1 - altFreq;
      }

      lookup[variant] = {
        variant,
        chrom,
        pos,
        pvalue,
        score,
        altAllele: alt,
        effectAllele: effect,
        altFreq: altFreq,
        effectFreq: effectFreq
      };
    });
    return lookup;
  }

  isAltEffect(variant_names) {  // Some calculations are sensitive to whether alt is the minor (effect) allele
    return variant_names.map(name => {
      const variant_data = this._variant_lookup[name];
      return variant_data.altAllele === variant_data.effectAllele;
    });
  }

  getEffectFreq(variant_names) {
    // Get the allele freq for the minor (effect) allele
    return variant_names.map(name => this._variant_lookup[name].effectFreq);
  }

  getScores(variant_names) {
    // Get single-variant scores
    return variant_names.map(name => this._variant_lookup[name].score);
  }

  getGroupVariants(variant_names) {
    // Return all that is known about a given set of variants
    return variant_names.map(name => this._variant_lookup[name]);
  }
}

// Utility class. Provides helper methods to access information about groups and generate subsets
class PortalGroupHelper {
  constructor(groups) {
    this._groups = groups;
    this._lookup = this._generateLookup(groups);
  }

  get data() {  // Raw unparsed data
    return this._groups;
  }

  byMask(selection) {  // str or array
    // Get all groups that identify as a specific category of mask- "limit the analysis to loss of function variants
    // in any gene"
    if (!Array.isArray(selection)) {
      selection = [selection]
    }
    selection = new Set(selection);

    const subset = this._groups.filter(group => selection.has(group.mask));
    return new this.constructor(subset);
  }

  byGroup(selection) {  // str or array
    // Get all groups based on a specific group name, regardless of mask. Eg, "all the ways to analyze data for a
    // given gene".
    if (!Array.isArray(selection)) {
      selection = [selection]
    }
    selection = new Set(selection);

    const subset = this._groups.filter(group => selection.has(group.group));
    return new this.constructor(subset);
  }

  _generateLookup(groups) {
    // We don't transform data, so this is a simple name -> position mapping
    return groups.reduce((acc, item, idx) => {
      const key = this._getKey(item.mask, item.group);
      acc[key] = idx;
      return acc;
    }, {});
  }

  _getKey(mask_name, group_name) {
    return `${mask_name},${group_name}`;
  }

  getOne(mask_name, group_name) {
    // Get a single group that is fully and uniquely identified by group + mask
    const key = this._getKey(mask_name, group_name);
    const pos = this._lookup[key];
    return this._groups[pos];
  }

  makeCovarianceMatrix(group, is_alt_effect) {
    // Helper method that expands the portal covariance format into a full matrix.
    // Load the covariance matrix from the response JSON
    const n_variants = group.variants.length;
    let covmat = new Array(n_variants);
    for (let i = 0; i < n_variants; i++) {
      covmat[i] = new Array(n_variants).fill(null);
    }

    let c = 0;
    for (let i = 0; i < n_variants; i++) {
      for (let j = i; j < n_variants; j++) {
        let v = group.covariance[c];
        let iAlt = is_alt_effect[i];
        let jAlt = is_alt_effect[j];

        /**
         * The API spec codes variant genotypes towards the alt allele. If the alt allele frequency
         * is > 0.5, that means we're not counting towards the minor (rare) allele, and we need to flip it around.
         * We don't flip when i == j because that element represents the variance of the variant's score, which is
         * invariant to which allele we code towards (but covariance is not.)
         *
         * We also don't flip when both the i variant and j variant need to be flipped (the ^ is XOR) because it would
         * just cancel out.
         */
        if (i !== j) {
          if ((!iAlt) ^ (!jAlt)) {
            v = -v;
          }
        }

        covmat[i][j] = v;
        covmat[j][i] = v;

        c += 1;
      }
    }

    covmat = numeric.mul(group.nSamples, covmat);
    return covmat;
  }
}

/**
 * Run one or more burden tests. This will operate in sequence: all specified tests on all specified masks
 *
 * The actual call signature of a burden test is pretty low-level. In addition to running the list of tests,
 *  this helper also restructures human-friendly mask and variant representations into a shape that works directly
 *  with the calculation.
 */
class PortalTestRunner {
  /**
   * Create a test runner object, using group and variant data of the form provided by `parsePortalJSON`. Generally,
   *  this helper is a convenience wrapper based on the raremetal.js API format spec, and hence it expects
   *  variant and group definitions to follow that spec.
   * @param groups PortalGroupHelper
   * @param variants PortalVariantsHelper
   * @param test_names {String[]|_AggregationTest[]}
   */
  constructor(groups, variants, test_names = []) {
    this.groups = groups;
    this.variants = variants;
    this._tests = [];

    test_names.forEach(name => this.addTest(name));
  }

  /**
   *
   * @param test {String|_AggregationTest}
   * @return {_AggregationTest}
   */
  addTest(test) {
    // Add a new test by name, or directly from an instance
    if (typeof test === 'string') {
      let type = AGGREGATION_TESTS[test];
      if (!type) {
        throw new Error(`Cannot make unknown test type: ${test}`);
      }
      test = new type.constructor();
    } else if (!(test instanceof _AggregationTest)) {
      throw new Error('Must specify test as name or instance');
    }
    this._tests.push(test);
    return test;
  }

  /**
   * Run every test on every group in the container and return results
   * @returns Promise A promise representing the fulfillment state of all tests being run
   */
  run() {
    let partials = [];

    this._tests.forEach(test => {
      this.groups.data.forEach(group => {
        partials.push(this._runOne.bind(this, test, group));
      });
    });
    // Despite the async syntax, ensure that each tests is run in series, to mitigate memory allocation errors when
    //  running many tests
    return partials.reduce((results, one_test) => {
      return results.then((all_prior) => {
        return one_test().then(one_res => {
          return [...all_prior, one_res];
        });
      });
    }, Promise.resolve([]));
  }

  /**
   *
   * @param {AggregationTest} test Instance for a single unit test
   * @param group {Object} Data corresponding to a specific group, following API format docs
   * @returns {{groupType: *, stat: *, test: *, pvalue: *, variants: (*|Array|string[]|Map), group: *, mask: (*|string|SVGMaskElement|string)}}
   * @private
   */
  _runOne(test, group) {
    // Helper method that translates portal data into the format expected by a test.
    const variants = group.variants;
    const scores = this.variants.getScores(variants);

    // Most calculations will require adjusting API data to ensure that minor allele is the effect allele
    const isAltEffect = this.variants.isAltEffect(variants);

    const cov = this.groups.makeCovarianceMatrix(group, isAltEffect);
    const mafs = this.variants.getEffectFreq(variants);
    let weights;  // TODO: The runner never actually uses the weights argument. Should it allow this?

    // Some test classes may return a raw value and others will return a promise. Wrap the result for consistency.
    let result = test.run(scores, cov, weights, mafs);
    return Promise.resolve(result)
      .then(([stat, pvalue]) => {
        // The results describe the group + several new fields for calculation results.
        return {
          groupType: group.groupType,
          group: group.group,
          mask: group.mask,
          variants: group.variants,

          test: test.key,
          stat,
          pvalue
        };
      });
  }

  /**
   * Generate a JSON representation of the results. Returns a Promise, because some methods may run asynchronously
   *  (eg via web workers), or require loading external libraries (eg webassembly)
   * @param results Array
   * @returns {Promise<{data: {groups: Promise<any> | Array, variants: *}} | never>}
   */
  toJSON(results) {
    // Output calculation results in a format that matches the "precomputed results" endpoint
    // By passing in an argument, user can format any set of results (even combining multiple runs)
    if (!results) {
      results = this.run();
    } else {
      results = Promise.resolve(results);
    }

    return results.then(group_results => {
      return {
        data: {
          variants: this.variants.data,
          groups: group_results,
        }
      }
    });
  }
}


/**
 * Define the runner code used by a web worker to run tests
 * Each inbound message initiates the process of running a test; each outbound message sends the results
 * @private
 */
function _make_runner(variants, test_names = []) {
  // The worker reuses the "run one" functionality of the basic test runner: one group/test pair at once
  // To level the load between workers, it doesn't know ahead of time what tasks will be assigned- just whatever is
  //  ready.

  const runner = new PortalTestRunner(null, variants, []);
  const tests = test_names.reduce((acc, name) => {
    acc[name] = runner.addTest(name);
  }, {});

  onmessage = function (event) { // Receive a test, execute, and notify when complete
    const [test_name, group] = event.data;
    runner._runOne(tests[test_name], group).then(result => {
      postMessage({ name: 'success', payload: result });
    });
  }
}

/**
 * A wrapper that allows parallel execution of tests using WebWorkers
 */
class PortalTestRunnerParallel {
  constructor(groups, variants, test_names = [], max_workers = 2) {
    if (typeof Worker === undefined) {
      throw new Error('Parallel test execution requires WebWorkers, which your environment does not support');
    }

    this.groups = groups;
    this.variants = variants;
    this._test_names = test_names;

    // Many tests are really fast to run, and a worker has some overhead.
    // Intelligently choose number of workers: only spin up separate threads if they are justified; eg more than
    //  one test is available to allocate each worker. We'll always run at least one worker, so a slow test wouldn't
    //  hang the browser.
    const num_tests = test_names.length * groups.data.length;
    this._num_workers = Math.min(max_workers, Math.ceil(num_tests / max_workers));
  }

  run() {
    // 1. Initiate the number of web workers requested (or run the executor directly and aggregate results)
    // 2. Post a message to each one describing which test we are running (proxy by name, b/c we can't pass objects into a worker)
    // 3. Receive the result of a test, and aggregate the results
    // 4. Close down each worker when it finishes its tasks and no new ones are available....
    // 5. Resolve the promise when the last task finishes (what's our criterion here? results.length === tests.length?)

    // Special case: if any single test throws an error, all tests are considered failed with an exception
    // The total run is computed from n tests * m groups
    const worker_code = new Blob(['(' + _make_runner.toString() + ')()'], { type: 'application/javascript' });
    const code_as_url = URL.createObjectURL(worker_code);
    let workers = [];

    return new Promise((resolve, reject) => {

      const all_tests = this._test_names.map(test => {
        this.groups.data.map(group => [test, group]);
      });
      const all_results = [];

      const _getNextTest = (worker) => {  // Send a new test to the web worker (if any are left)
        const spec = all_tests.pop();
        if (spec) {
          worker.postMessage(spec);
        }
      };

      // Generate n independent workers, each one equipped to start tasks
      workers = Array.from(new Array(this._num_workers), () => {
        const worker = new Worker(code_as_url);
        worker.onmessage = event => { // Each time a result is received, give the worker a new task
          all_results.push(event.data);

          if (all_results.length === all_tests.length) {
            // Assumption: if what we are receiving is the last expected result, it's ok to end calc and stop all
            //  web workers, because they already finished their calculations
            resolve(all_results);
          } else {
            _getNextTest(worker);
          }
        };
        // If any single worker fails, we consider the entire calculation failed and end the runner with this error
        worker.onerror = e => reject(e);
        return worker;
      });
      workers.forEach(worker => _getNextTest(worker));
    }).finally(() => {
      // Cleanup to prevent memory leaks: terminate all web workers and revoke code url to prevent memory leaks
      URL.revokeObjectURL(code_as_url);
      workers.forEach(worker => worker.terminate());
    });
  }
}


function parsePortalJSON(json) {
  const data = json.data || json;
  const groups = new PortalGroupHelper(data.groups.map(item => {
    // Each group should have access to fields that, in portal json, are defined once globally
    item.nSamples = data.nSamples;
    item.sigmaSquared = data.sigmaSquared;
    return item;
  }));
  const variants = new PortalVariantsHelper(data.variants);
  return [groups, variants];
}

export { PortalVariantsHelper as _PortalVariantsHelper, PortalGroupHelper as _PortalGroupHelper }; // testing only

export { parsePortalJSON, PortalTestRunner, PortalTestRunnerParallel, AGGREGATION_TESTS };
