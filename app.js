/* ============================================================
   NC Interpolation — Interpolation Calculator
   Lagrange, Newton Divided/Forward/Backward, Stirling
   ============================================================ */

(function () {
  'use strict';

  // ===== FRACTION ARITHMETIC =====
  function gcd(a, b) {
    a = Math.abs(a); b = Math.abs(b);
    while (b) { [a, b] = [b, a % b]; }
    return a || 1;
  }

  class Frac {
    constructor(num, den = 1) {
      if (den === 0) throw new Error('Division by zero');
      if (den < 0) { num = -num; den = -den; }
      const g = gcd(Math.abs(num), Math.abs(den));
      this.n = num / g;
      this.d = den / g;
    }
    add(o) { return new Frac(this.n * o.d + o.n * this.d, this.d * o.d); }
    sub(o) { return new Frac(this.n * o.d - o.n * this.d, this.d * o.d); }
    mul(o) { return new Frac(this.n * o.n, this.d * o.d); }
    div(o) {
      if (o.n === 0) throw new Error('Division by zero');
      return new Frac(this.n * o.d, this.d * o.n);
    }
    neg() { return new Frac(-this.n, this.d); }
    isZero() { return this.n === 0; }
    toNumber() { return this.n / this.d; }
    toString() { return this.d === 1 ? '' + this.n : this.n + '/' + this.d; }
    toHTML() {
      if (this.d === 1) return '<span class="num-int">' + this.n + '</span>';
      var sign = this.n < 0 ? '<span class="num-int">\u2212</span>' : '';
      return sign + '<span class="frac"><span class="frac-num">' + Math.abs(this.n) + '</span><span class="frac-den">' + this.d + '</span></span>';
    }
    static zero() { return new Frac(0); }
    static one() { return new Frac(1); }
    static from(val) {
      if (val instanceof Frac) return new Frac(val.n, val.d);
      if (Number.isInteger(val)) return new Frac(val);
      return Frac._fromFloat(val);
    }
    static _fromFloat(x) {
      if (Math.abs(x - Math.round(x)) < 1e-9) return new Frac(Math.round(x));
      var neg = x < 0; x = Math.abs(x);
      var h1 = 1, h2 = 0, k1 = 0, k2 = 1, b = x;
      for (var i = 0; i < 25; i++) {
        var a = Math.floor(b);
        var h = a * h1 + h2, k = a * k1 + k2;
        if (k > 100000) break;
        h2 = h1; h1 = h; k2 = k1; k1 = k;
        if (Math.abs(b - a) < 1e-10) break;
        b = 1 / (b - a);
      }
      return new Frac(neg ? -h1 : h1, k1);
    }
  }

  // ===== HELPERS =====
  function fmt(v) {
    if (typeof v === 'number') {
      if (Number.isInteger(v)) return v.toString();
      return parseFloat(v.toFixed(8)).toString();
    }
    return v.toString();
  }

  function fmtHTML(v) {
    if (v instanceof Frac) return v.toHTML();
    if (typeof v === 'number') {
      if (Number.isInteger(v)) return '<span class="num-int">' + v + '</span>';
      return '<span class="num-int">' + parseFloat(v.toFixed(8)) + '</span>';
    }
    return String(v);
  }

  function escapeHTML(s) { return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;'); }

  // Check if data is equally spaced
  function isEquallySpaced(xs) {
    if (xs.length < 2) return { equal: true, h: 0 };
    var h = xs[1] - xs[0];
    for (var i = 2; i < xs.length; i++) {
      if (Math.abs((xs[i] - xs[i-1]) - h) > 1e-9) return { equal: false, h: 0 };
    }
    return { equal: true, h: h };
  }

  // ===== DIFFERENCE TABLES =====

  // Divided difference table: table[i][j] = f[x_i,...,x_{i+j}]
  function dividedDiffTable(xs, ys) {
    var n = xs.length;
    var table = [];
    for (var i = 0; i < n; i++) {
      table[i] = [ys[i]];
    }
    for (var j = 1; j < n; j++) {
      for (var i = 0; i < n - j; i++) {
        table[i][j] = (table[i+1][j-1] - table[i][j-1]) / (xs[i+j] - xs[i]);
      }
    }
    return table;
  }

  // Forward difference table: delta[i][j] = Δ^j f(x_i)
  function forwardDiffTable(ys) {
    var n = ys.length;
    var table = [];
    for (var i = 0; i < n; i++) table[i] = [ys[i]];
    for (var j = 1; j < n; j++) {
      for (var i = 0; i < n - j; i++) {
        table[i][j] = table[i+1][j-1] - table[i][j-1];
      }
    }
    return table;
  }

  // ===== INTERPOLATION METHODS =====

  // 1. Lagrange Interpolation
  function lagrangeInterpolation(xs, ys, xVal) {
    var n = xs.length;
    var result = 0;
    var steps = [];
    var terms = [];

    for (var i = 0; i < n; i++) {
      var num = 1, den = 1;
      var numParts = [], denParts = [];
      for (var j = 0; j < n; j++) {
        if (i !== j) {
          num *= (xVal - xs[j]);
          den *= (xs[i] - xs[j]);
          numParts.push('(' + fmt(xVal) + ' - ' + fmt(xs[j]) + ')');
          denParts.push('(' + fmt(xs[i]) + ' - ' + fmt(xs[j]) + ')');
        }
      }
      var Li = num / den;
      var term = Li * ys[i];
      result += term;
      terms.push({ i: i, Li: Li, yi: ys[i], term: term, numParts: numParts, denParts: denParts });
      steps.push(
        'L' + i + '(x) = ' + numParts.join(' × ') + ' / ' + denParts.join(' × ') +
        ' = ' + fmt(Li)
      );
      steps.push('  L' + i + ' × f(x' + i + ') = ' + fmt(Li) + ' × ' + fmt(ys[i]) + ' = ' + fmt(term));
    }
    steps.push('');
    steps.push('f(' + fmt(xVal) + ') = ' + terms.map(function(t){ return fmt(t.term); }).join(' + ') + ' = ' + fmt(result));

    return { value: result, steps: steps, terms: terms };
  }

  // Build Lagrange polynomial symbolically (simplified coefficients)
  function lagrangePolynomial(xs, ys) {
    var n = xs.length;
    // Represent polynomial as array of coefficients [c0, c1, ..., cn-1]
    var poly = new Array(n).fill(0);
    for (var i = 0; i < n; i++) {
      var den = 1;
      for (var j = 0; j < n; j++) {
        if (i !== j) den *= (xs[i] - xs[j]);
      }
      // Build product polynomial (x - x0)(x - x1)... excluding xi
      var basis = [1]; // start with 1
      for (var j = 0; j < n; j++) {
        if (i !== j) {
          var newBasis = new Array(basis.length + 1).fill(0);
          for (var k = 0; k < basis.length; k++) {
            newBasis[k + 1] += basis[k];        // x term
            newBasis[k] += basis[k] * (-xs[j]);  // constant term
          }
          basis = newBasis;
        }
      }
      var coeff = ys[i] / den;
      for (var k = 0; k < basis.length; k++) {
        poly[k] += coeff * basis[k];
      }
    }
    return poly;
  }

  // 2. Newton's Divided Difference
  function newtonDividedDiff(xs, ys, xVal) {
    var n = xs.length;
    var table = dividedDiffTable(xs, ys);
    var result = table[0][0];
    var product = 1;
    var steps = [];
    var coeffs = [table[0][0]];

    steps.push('Divided difference table coefficients (top diagonal):');
    for (var j = 0; j < n; j++) {
      steps.push('  f[' + xs.slice(0, j+1).map(fmt).join(', ') + '] = ' + fmt(table[0][j]));
    }
    steps.push('');
    steps.push('P(x) = f[x₀] + f[x₀,x₁](x−x₀) + f[x₀,x₁,x₂](x−x₀)(x−x₁) + ...');
    steps.push('');

    var termStrs = [fmt(table[0][0])];
    for (var j = 1; j < n; j++) {
      product *= (xVal - xs[j-1]);
      var term = table[0][j] * product;
      result += term;
      coeffs.push(table[0][j]);

      var factors = [];
      for (var k = 0; k < j; k++) factors.push('(' + fmt(xVal) + ' - ' + fmt(xs[k]) + ')');
      steps.push('Term ' + j + ': ' + fmt(table[0][j]) + ' × ' + factors.join(' × ') + ' = ' + fmt(term));
      termStrs.push(fmt(term));
    }
    steps.push('');
    steps.push('f(' + fmt(xVal) + ') = ' + termStrs.join(' + ') + ' = ' + fmt(result));

    return { value: result, steps: steps, table: table, coeffs: coeffs };
  }

  // 3. Newton's Forward Difference
  function newtonForwardDiff(xs, ys, xVal) {
    var n = xs.length;
    var h = xs[1] - xs[0];
    var table = forwardDiffTable(ys);
    var s = (xVal - xs[0]) / h;
    var result = table[0][0];
    var steps = [];

    steps.push('Step size h = ' + fmt(h));
    steps.push('s = (x − x₀) / h = (' + fmt(xVal) + ' − ' + fmt(xs[0]) + ') / ' + fmt(h) + ' = ' + fmt(s));
    steps.push('');
    steps.push('Forward differences (Δ⁰f₀, Δ¹f₀, Δ²f₀, ...):');
    for (var j = 0; j < n; j++) {
      steps.push('  Δ' + (j === 0 ? '⁰' : superscript(j)) + 'f₀ = ' + fmt(table[0][j]));
    }
    steps.push('');

    var sTerm = 1;
    var factorial = 1;
    var termStrs = [fmt(table[0][0])];

    for (var j = 1; j < n; j++) {
      sTerm *= (s - (j - 1));
      factorial *= j;
      var term = (sTerm / factorial) * table[0][j];
      result += term;
      steps.push('Term ' + j + ': s' + productTermDesc(j, s) + ' / ' + j + '! × Δ' + superscript(j) + 'f₀ = ' + fmt(sTerm/factorial) + ' × ' + fmt(table[0][j]) + ' = ' + fmt(term));
      termStrs.push(fmt(term));
    }
    steps.push('');
    steps.push('f(' + fmt(xVal) + ') = ' + termStrs.join(' + ') + ' = ' + fmt(result));

    return { value: result, steps: steps, table: table, h: h, s: s };
  }

  // 4. Newton's Backward Difference
  function newtonBackwardDiff(xs, ys, xVal) {
    var n = xs.length;
    var h = xs[1] - xs[0];
    // Backward diff table: ∇^j f(x_n) uses the LAST value
    var table = forwardDiffTable(ys); // same table, different reading
    var xn = xs[n - 1];
    var s = (xVal - xn) / h;
    var result = ys[n - 1];
    var steps = [];

    steps.push('Step size h = ' + fmt(h));
    steps.push('s = (x − xₙ) / h = (' + fmt(xVal) + ' − ' + fmt(xn) + ') / ' + fmt(h) + ' = ' + fmt(s));
    steps.push('');
    steps.push('Backward differences (∇⁰fₙ, ∇¹fₙ, ∇²fₙ, ...):');

    // Backward differences: ∇^j f_n = table[n-1-j][j] (diagonal from bottom-left)
    var backDiffs = [ys[n-1]];
    for (var j = 1; j < n; j++) {
      backDiffs.push(table[n - 1 - j][j]);
      steps.push('  ∇' + superscript(j) + 'fₙ = ' + fmt(table[n - 1 - j][j]));
    }
    steps.push('');

    var sTerm = 1;
    var factorial = 1;
    var termStrs = [fmt(ys[n-1])];

    for (var j = 1; j < n; j++) {
      sTerm *= (s + (j - 1));
      factorial *= j;
      var term = (sTerm / factorial) * backDiffs[j];
      result += term;
      steps.push('Term ' + j + ': s' + backProductTermDesc(j, s) + ' / ' + j + '! × ∇' + superscript(j) + 'fₙ = ' + fmt(sTerm/factorial) + ' × ' + fmt(backDiffs[j]) + ' = ' + fmt(term));
      termStrs.push(fmt(term));
    }
    steps.push('');
    steps.push('f(' + fmt(xVal) + ') = ' + termStrs.join(' + ') + ' = ' + fmt(result));

    return { value: result, steps: steps, table: table, h: h, s: s, backDiffs: backDiffs };
  }

  // 5. Stirling's Formula
  function stirlingInterpolation(xs, ys, xVal) {
    var n = xs.length;
    if (n % 2 === 0) throw new Error("Stirling's formula requires an odd number of data points (central point needed).");
    var h = xs[1] - xs[0];
    var table = forwardDiffTable(ys);
    var mid = Math.floor(n / 2);
    var x0 = xs[mid];
    var s = (xVal - x0) / h;
    var result = ys[mid];
    var steps = [];

    steps.push('Central point x₀ = ' + fmt(x0) + ' (index ' + mid + ')');
    steps.push('Step size h = ' + fmt(h));
    steps.push('s = (x − x₀) / h = (' + fmt(xVal) + ' − ' + fmt(x0) + ') / ' + fmt(h) + ' = ' + fmt(s));
    steps.push('');

    // Stirling uses averages of forward differences at different levels
    // δ^(2k) f₀ = Δ^(2k) f_{mid-k}  (even central differences)
    // For odd: mean of Δ^(2k+1) f_{mid-k-1} and Δ^(2k+1) f_{mid-k}

    var termStrs = [fmt(ys[mid])];

    // We'll compute up to available order
    var maxOrder = n - 1;
    var sTerm_odd = s;  // s for first odd term
    var sTerm_even = s * s; // s² for first even term
    var factorial_odd = 1;
    var factorial_even = 2;

    for (var k = 1; k <= maxOrder; k++) {
      if (k % 2 === 1) {
        // Odd order: use mean of two adjacent differences
        var order = k;
        var idx1 = mid - Math.floor(order / 2) - 1;
        var idx2 = mid - Math.floor(order / 2);
        if (idx1 < 0 || idx2 < 0 || !table[idx1] || !table[idx2] ||
            table[idx1][order] === undefined || table[idx2][order] === undefined) break;

        var mean = (table[idx1][order] + table[idx2][order]) / 2;
        // s product: s for k=1, s(s²-1) for k=3, s(s²-1)(s²-4) for k=5, etc.
        var halfK = Math.floor(k / 2);
        if (k === 1) {
          sTerm_odd = s;
          factorial_odd = 1;
        } else {
          // k=3: multiply by (s²-1)/( 2*3 )
          sTerm_odd *= (s * s - halfK * halfK + (halfK-1)*(halfK-1) === undefined ? 1 : 1);
          // Re-derive: for Stirling, odd term k uses: s * prod_{j=1}^{(k-1)/2} (s²-j²) / k!
        }
        // Simpler: compute directly
        var prod = s;
        var fact = 1;
        for (var j = 1; j <= halfK; j++) {
          prod *= (s * s - j * j);
        }
        for (var j = 1; j <= k; j++) fact *= j;

        var term = (prod / fact) * mean;
        result += term;
        steps.push('Odd term (order ' + k + '): μδ' + superscript(k) + 'f₀ = (' + fmt(table[idx1][order]) + ' + ' + fmt(table[idx2][order]) + ')/2 = ' + fmt(mean));
        steps.push('  Coefficient = ' + fmt(prod/fact) + ', term = ' + fmt(term));
        termStrs.push(fmt(term));
      } else {
        // Even order: central difference
        var order = k;
        var idx = mid - k / 2;
        if (idx < 0 || !table[idx] || table[idx][order] === undefined) break;

        var diff = table[idx][order];
        var halfK = k / 2;
        var prod = s * s;
        var fact = 2;
        for (var j = 1; j < halfK; j++) {
          prod *= (s * s - j * j);
        }
        for (var j = 3; j <= k; j++) fact *= j;

        var term = (prod / fact) * diff;
        result += term;
        steps.push('Even term (order ' + k + '): δ' + superscript(k) + 'f₀ = ' + fmt(diff));
        steps.push('  Coefficient = ' + fmt(prod/fact) + ', term = ' + fmt(term));
        termStrs.push(fmt(term));
      }
    }

    steps.push('');
    steps.push('f(' + fmt(xVal) + ') = ' + termStrs.join(' + ') + ' = ' + fmt(result));

    return { value: result, steps: steps, table: table, h: h, s: s, mid: mid };
  }

  // ===== HELPER: superscript numbers =====
  function superscript(n) {
    var sup = '⁰¹²³⁴⁵⁶⁷⁸⁹';
    return String(n).split('').map(function(c) { return sup[parseInt(c)] || c; }).join('');
  }

  function productTermDesc(j, s) {
    var parts = [];
    for (var k = 0; k < j; k++) {
      if (k === 0) parts.push('');
      else parts.push('(s − ' + k + ')');
    }
    return parts.filter(Boolean).length ? ' × ' + parts.filter(Boolean).join(' × ') : '';
  }

  function backProductTermDesc(j, s) {
    var parts = [];
    for (var k = 0; k < j; k++) {
      if (k === 0) parts.push('');
      else parts.push('(s + ' + k + ')');
    }
    return parts.filter(Boolean).length ? ' × ' + parts.filter(Boolean).join(' × ') : '';
  }

  // Format polynomial from coefficients array
  function formatPolynomial(coeffs) {
    var terms = [];
    for (var i = coeffs.length - 1; i >= 0; i--) {
      var c = coeffs[i];
      if (Math.abs(c) < 1e-10) continue;
      var cs = '';
      if (i === 0) {
        cs = fmt(c);
      } else if (Math.abs(c - 1) < 1e-10) {
        cs = '';
      } else if (Math.abs(c + 1) < 1e-10) {
        cs = '-';
      } else {
        cs = fmt(c);
      }
      var xpart = '';
      if (i === 1) xpart = 'x';
      else if (i > 1) xpart = 'x' + superscript(i);

      var termStr = cs + xpart;
      if (!termStr || termStr === '-') termStr = cs + xpart || fmt(c);
      terms.push(termStr);
    }
    if (terms.length === 0) return '0';
    var result = terms[0];
    for (var i = 1; i < terms.length; i++) {
      if (terms[i].startsWith('-')) {
        result += ' − ' + terms[i].substring(1);
      } else {
        result += ' + ' + terms[i];
      }
    }
    return 'P(x) = ' + result;
  }

  // ===== LaTeX Generation =====
  function generateLaTeX(method, xs, ys, xVal, resultObj) {
    var lines = [];
    lines.push('% Interpolation — ' + methodName(method));
    lines.push('\\text{Given data points:}');
    lines.push('\\begin{array}{c|c}');
    lines.push('x & f(x) \\\\ \\hline');
    for (var i = 0; i < xs.length; i++) {
      lines.push(fmt(xs[i]) + ' & ' + fmt(ys[i]) + ' \\\\');
    }
    lines.push('\\end{array}');
    lines.push('');
    lines.push('\\text{Find } f(' + fmt(xVal) + ')');
    lines.push('');
    lines.push('f(' + fmt(xVal) + ') = ' + fmt(resultObj.value));
    return lines.join('\n');
  }

  // ===== METHOD NAMES =====
  function methodName(m) {
    var names = {
      'lagrange': 'Lagrange Interpolation',
      'newton-divided': "Newton's Divided Difference",
      'newton-forward': "Newton's Forward Difference",
      'newton-backward': "Newton's Backward Difference",
      'stirling': "Stirling's Formula",
      'auto': 'Auto-Detect'
    };
    return names[m] || m;
  }

  // ===== PLOTTING =====
  function plotInterpolation(canvas, xs, ys, xVal, yVal, polyEval) {
    var ctx = canvas.getContext('2d');
    var W = canvas.width = canvas.offsetWidth * 2;
    var H = canvas.height = 400 * 2;
    ctx.scale(2, 2);
    var w = W / 2, h = H / 2;
    var isDark = document.documentElement.getAttribute('data-theme') === 'dark';

    // Clear
    ctx.fillStyle = isDark ? '#161830' : '#f4f6f9';
    ctx.fillRect(0, 0, w, h);

    // Compute plot bounds
    var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
    var yMin = Math.min.apply(null, ys), yMax = Math.max.apply(null, ys);
    if (xVal !== null) {
      xMin = Math.min(xMin, xVal);
      xMax = Math.max(xMax, xVal);
    }
    var xPad = (xMax - xMin) * 0.15 || 1;
    var yPad = (yMax - yMin) * 0.15 || 1;
    xMin -= xPad; xMax += xPad;
    yMin -= yPad; yMax += yPad;

    // Evaluate curve at many points for yMin/yMax
    var curvePoints = [];
    var steps = 200;
    for (var i = 0; i <= steps; i++) {
      var cx = xMin + (xMax - xMin) * i / steps;
      var cy = polyEval(cx);
      if (isFinite(cy)) {
        curvePoints.push([cx, cy]);
        yMin = Math.min(yMin, cy);
        yMax = Math.max(yMax, cy);
      }
    }
    yPad = (yMax - yMin) * 0.1 || 1;
    yMin -= yPad; yMax += yPad;

    var pad = { left: 50, right: 20, top: 20, bottom: 35 };
    var pw = w - pad.left - pad.right;
    var ph = h - pad.top - pad.bottom;

    function toScreen(px, py) {
      return [
        pad.left + (px - xMin) / (xMax - xMin) * pw,
        pad.top + (1 - (py - yMin) / (yMax - yMin)) * ph
      ];
    }

    // Grid
    ctx.strokeStyle = isDark ? '#2e3158' : '#e2e4ea';
    ctx.lineWidth = 0.5;
    ctx.font = '10px -apple-system, sans-serif';
    ctx.fillStyle = isDark ? '#6b6f90' : '#8b8ea8';
    ctx.textAlign = 'center';

    // X axis labels
    var xTicks = 6;
    for (var i = 0; i <= xTicks; i++) {
      var tx = xMin + (xMax - xMin) * i / xTicks;
      var sp = toScreen(tx, yMin);
      ctx.beginPath();
      ctx.moveTo(sp[0], pad.top);
      ctx.lineTo(sp[0], pad.top + ph);
      ctx.stroke();
      ctx.fillText(tx.toFixed(2), sp[0], h - 5);
    }

    // Y axis labels
    ctx.textAlign = 'right';
    var yTicks = 5;
    for (var i = 0; i <= yTicks; i++) {
      var ty = yMin + (yMax - yMin) * i / yTicks;
      var sp = toScreen(xMin, ty);
      ctx.beginPath();
      ctx.moveTo(pad.left, sp[1]);
      ctx.lineTo(pad.left + pw, sp[1]);
      ctx.stroke();
      ctx.fillText(ty.toFixed(2), pad.left - 5, sp[1] + 3);
    }

    // Axes
    ctx.strokeStyle = isDark ? '#a0a3c0' : '#5a5d7a';
    ctx.lineWidth = 1;
    // X axis (y=0 if in range)
    if (yMin <= 0 && yMax >= 0) {
      var ay = toScreen(0, 0);
      ctx.beginPath();
      ctx.moveTo(pad.left, ay[1]);
      ctx.lineTo(pad.left + pw, ay[1]);
      ctx.stroke();
    }
    // Y axis (x=0 if in range)
    if (xMin <= 0 && xMax >= 0) {
      var ax = toScreen(0, 0);
      ctx.beginPath();
      ctx.moveTo(ax[0], pad.top);
      ctx.lineTo(ax[0], pad.top + ph);
      ctx.stroke();
    }

    // Curve
    ctx.strokeStyle = isDark ? '#4cc9f0' : '#4361ee';
    ctx.lineWidth = 2;
    ctx.beginPath();
    var started = false;
    for (var i = 0; i < curvePoints.length; i++) {
      var sp = toScreen(curvePoints[i][0], curvePoints[i][1]);
      if (sp[0] >= pad.left && sp[0] <= pad.left + pw && sp[1] >= pad.top && sp[1] <= pad.top + ph) {
        if (!started) { ctx.moveTo(sp[0], sp[1]); started = true; }
        else ctx.lineTo(sp[0], sp[1]);
      } else {
        started = false;
      }
    }
    ctx.stroke();

    // Data points
    ctx.fillStyle = isDark ? '#4cc9f0' : '#4361ee';
    for (var i = 0; i < xs.length; i++) {
      var sp = toScreen(xs[i], ys[i]);
      ctx.beginPath();
      ctx.arc(sp[0], sp[1], 5, 0, Math.PI * 2);
      ctx.fill();
      ctx.strokeStyle = isDark ? '#1c1f3a' : '#ffffff';
      ctx.lineWidth = 2;
      ctx.stroke();
    }

    // Interpolated point
    if (xVal !== null && yVal !== null && isFinite(yVal)) {
      var sp = toScreen(xVal, yVal);
      ctx.fillStyle = '#e74c3c';
      ctx.beginPath();
      ctx.arc(sp[0], sp[1], 7, 0, Math.PI * 2);
      ctx.fill();
      ctx.strokeStyle = isDark ? '#1c1f3a' : '#ffffff';
      ctx.lineWidth = 2;
      ctx.stroke();

      ctx.fillStyle = isDark ? '#e4e6f0' : '#1a1b2e';
      ctx.font = 'bold 11px -apple-system, sans-serif';
      ctx.textAlign = 'left';
      ctx.fillText('(' + fmt(xVal) + ', ' + parseFloat(yVal.toFixed(6)) + ')', sp[0] + 10, sp[1] - 8);
    }
  }

  // ===== BUILD DIFFERENCE TABLE HTML =====
  function buildDiffTableHTML(xs, ys, table, type) {
    var n = xs.length;
    var html = '<table class="diff-table"><thead><tr>';
    html += '<th>x</th><th>f(x)</th>';

    if (type === 'divided') {
      for (var j = 1; j < n; j++) {
        html += '<th>Order ' + j + '</th>';
      }
    } else {
      for (var j = 1; j < n; j++) {
        html += '<th>Δ' + superscript(j) + 'f</th>';
      }
    }
    html += '</tr></thead><tbody>';

    for (var i = 0; i < n; i++) {
      html += '<tr>';
      html += '<td>' + fmt(xs[i]) + '</td>';
      for (var j = 0; j < n; j++) {
        if (table[i] && table[i][j] !== undefined) {
          html += '<td>' + parseFloat(table[i][j].toFixed(8)) + '</td>';
        } else {
          html += '<td></td>';
        }
      }
      html += '</tr>';
    }
    html += '</tbody></table>';
    return html;
  }

  // ===== SHARE LINK =====
  function buildShareURL(xs, ys, xVal, method) {
    var params = new URLSearchParams();
    params.set('x', xs.join(','));
    params.set('y', ys.join(','));
    if (xVal !== null && xVal !== '') params.set('at', xVal);
    params.set('m', method);
    return location.origin + location.pathname + '?' + params.toString();
  }

  function loadFromURL() {
    var params = new URLSearchParams(location.search);
    if (!params.has('x') || !params.has('y')) return false;
    var xs = params.get('x').split(',').map(Number);
    var ys = params.get('y').split(',').map(Number);
    if (xs.length !== ys.length || xs.length < 2) return false;

    var numPts = document.getElementById('num-points');
    numPts.value = xs.length;
    buildDataTable(xs.length);

    var xInputs = document.querySelectorAll('.x-input');
    var yInputs = document.querySelectorAll('.y-input');
    for (var i = 0; i < xs.length; i++) {
      if (xInputs[i]) xInputs[i].value = xs[i];
      if (yInputs[i]) yInputs[i].value = ys[i];
    }

    if (params.has('at')) {
      document.getElementById('interp-x').value = params.get('at');
    }
    if (params.has('m')) {
      document.getElementById('method-select').value = params.get('m');
    }
    return true;
  }

  // ===== HISTORY =====
  var HISTORY_KEY = 'nc-interp-history';
  var MAX_HISTORY = 20;

  function loadHistory() {
    try { return JSON.parse(localStorage.getItem(HISTORY_KEY)) || []; }
    catch(e) { return []; }
  }

  function saveHistory(entry) {
    var hist = loadHistory();
    hist.unshift(entry);
    if (hist.length > MAX_HISTORY) hist = hist.slice(0, MAX_HISTORY);
    localStorage.setItem(HISTORY_KEY, JSON.stringify(hist));
    renderHistory();
  }

  function clearHistory() {
    localStorage.removeItem(HISTORY_KEY);
    renderHistory();
  }

  function renderHistory() {
    var list = document.getElementById('history-list');
    var hist = loadHistory();
    if (hist.length === 0) {
      list.innerHTML = '<p class="history-empty">No calculations yet</p>';
      return;
    }
    var html = '';
    hist.forEach(function(item, idx) {
      html += '<div class="history-item" data-idx="' + idx + '">';
      html += '<div class="history-meta">';
      html += '<span class="history-method">' + escapeHTML(methodName(item.method)) + '</span>';
      html += '<span class="history-size">' + item.n + ' pts</span>';
      html += '<span class="history-date">' + escapeHTML(item.date) + '</span>';
      html += '</div>';
      html += '<div class="history-preview">f(' + escapeHTML(fmt(item.xVal)) + ') = ' + escapeHTML(fmt(item.result)) + '</div>';
      html += '</div>';
    });
    list.innerHTML = html;

    list.querySelectorAll('.history-item').forEach(function(el) {
      el.addEventListener('click', function() {
        var idx = parseInt(el.getAttribute('data-idx'));
        var item = loadHistory()[idx];
        if (!item) return;
        restoreFromHistory(item);
      });
    });
  }

  function restoreFromHistory(item) {
    var numPts = document.getElementById('num-points');
    numPts.value = item.n;
    buildDataTable(item.n);
    var xInputs = document.querySelectorAll('.x-input');
    var yInputs = document.querySelectorAll('.y-input');
    for (var i = 0; i < item.xs.length; i++) {
      if (xInputs[i]) xInputs[i].value = item.xs[i];
      if (yInputs[i]) yInputs[i].value = item.ys[i];
    }
    document.getElementById('interp-x').value = item.xVal;
    document.getElementById('method-select').value = item.method;
    calculate();
  }

  // ===== BUILD DATA TABLE =====
  function buildDataTable(n) {
    var body = document.getElementById('data-table-body');
    body.innerHTML = '';
    for (var i = 0; i < n; i++) {
      var tr = document.createElement('tr');
      tr.innerHTML = '<td class="row-num">' + i + '</td>' +
        '<td><input type="text" class="data-input x-input" placeholder="x' + i + '" data-idx="' + i + '"></td>' +
        '<td><input type="text" class="data-input y-input" placeholder="f(x' + i + ')" data-idx="' + i + '"></td>';
      body.appendChild(tr);
    }
  }

  // ===== CALCULATE =====
  function calculate() {
    var output = document.getElementById('output');
    output.innerHTML = '';

    // Read inputs
    var xInputs = document.querySelectorAll('.x-input');
    var yInputs = document.querySelectorAll('.y-input');
    var xs = [], ys = [];
    var hasError = false;

    xInputs.forEach(function(inp) {
      inp.classList.remove('input-error');
      var v = parseFloat(inp.value);
      if (isNaN(v)) { inp.classList.add('input-error'); hasError = true; }
      xs.push(v);
    });
    yInputs.forEach(function(inp) {
      inp.classList.remove('input-error');
      var v = parseFloat(inp.value);
      if (isNaN(v)) { inp.classList.add('input-error'); hasError = true; }
      ys.push(v);
    });

    var xValInput = document.getElementById('interp-x');
    var xVal = parseFloat(xValInput.value);
    if (isNaN(xVal)) {
      xValInput.classList.add('input-error');
      hasError = true;
    } else {
      xValInput.classList.remove('input-error');
    }

    if (hasError) {
      output.innerHTML = '<div class="error-message">Please fill in all data points and the interpolation value with valid numbers.</div>';
      return;
    }

    // Check for duplicate x values
    var xSet = new Set(xs);
    if (xSet.size !== xs.length) {
      output.innerHTML = '<div class="error-message">Duplicate x-values detected. Each x must be unique.</div>';
      return;
    }

    // Sort by x
    var pairs = xs.map(function(x, i) { return [x, ys[i]]; });
    pairs.sort(function(a, b) { return a[0] - b[0]; });
    xs = pairs.map(function(p) { return p[0]; });
    ys = pairs.map(function(p) { return p[1]; });

    var method = document.getElementById('method-select').value;
    var spacing = isEquallySpaced(xs);
    var compareMode = document.getElementById('compare-toggle').checked;

    // Auto-detect method
    if (method === 'auto') {
      if (spacing.equal) {
        // Choose based on position of xVal
        if (xVal <= xs[0] + (xs[xs.length-1] - xs[0]) * 0.4) {
          method = 'newton-forward';
        } else if (xVal >= xs[0] + (xs[xs.length-1] - xs[0]) * 0.6) {
          method = 'newton-backward';
        } else if (xs.length % 2 === 1) {
          method = 'stirling';
        } else {
          method = 'newton-forward';
        }
      } else {
        method = 'newton-divided';
      }
    }

    // Validate method vs spacing
    var needsEqual = ['newton-forward', 'newton-backward', 'stirling'];
    if (needsEqual.indexOf(method) !== -1 && !spacing.equal) {
      output.innerHTML = '<div class="error-message">' + methodName(method) + ' requires equally spaced data points. Your data is not equally spaced. Try Lagrange or Newton\'s Divided Difference instead.</div>';
      return;
    }

    // Show spacing note
    var note = document.getElementById('spacing-note');
    if (spacing.equal) {
      note.style.display = 'block';
      note.textContent = 'Data is equally spaced (h = ' + fmt(spacing.h) + '). All methods available.';
    } else {
      note.style.display = 'block';
      note.textContent = 'Data is NOT equally spaced. Only Lagrange and Newton\'s Divided Difference are applicable.';
    }

    function runMethod(m) {
      try {
        switch (m) {
          case 'lagrange': return lagrangeInterpolation(xs, ys, xVal);
          case 'newton-divided': return newtonDividedDiff(xs, ys, xVal);
          case 'newton-forward': return newtonForwardDiff(xs, ys, xVal);
          case 'newton-backward': return newtonBackwardDiff(xs, ys, xVal);
          case 'stirling': return stirlingInterpolation(xs, ys, xVal);
          default: return lagrangeInterpolation(xs, ys, xVal);
        }
      } catch (e) {
        return { error: e.message };
      }
    }

    function buildResultHTML(m, res) {
      if (res.error) {
        return '<div class="error-message">' + escapeHTML(res.error) + '</div>';
      }

      var html = '';
      html += '<div class="result-panel">';
      html += '<h3 class="result-title">' + escapeHTML(methodName(m)) + '</h3>';

      // Result value
      html += '<div class="result-value">f(' + fmt(xVal) + ') = ' + parseFloat(res.value.toFixed(10)) + '</div>';

      // Difference table
      if (res.table) {
        html += '<div class="result-section"><h3>Difference Table</h3>';
        var tableType = (m === 'newton-divided') ? 'divided' : 'forward';
        html += buildDiffTableHTML(xs, ys, res.table, tableType);
        html += '</div>';
      }

      // Polynomial (for Lagrange)
      if (m === 'lagrange') {
        try {
          var poly = lagrangePolynomial(xs, ys);
          html += '<div class="result-section"><h3>Interpolating Polynomial</h3>';
          html += '<div class="polynomial-display">' + escapeHTML(formatPolynomial(poly)) + '</div>';
          html += '</div>';
        } catch(e) {}
      }

      // Steps
      if (res.steps && res.steps.length) {
        html += '<details class="steps-details"><summary>Step-by-Step Solution</summary>';
        html += '<div class="steps-content">';
        res.steps.forEach(function(s) {
          if (s === '') {
            html += '<br>';
          } else if (s.startsWith('  ')) {
            html += '<div class="step-block"><div class="step-line">' + escapeHTML(s.trim()) + '</div></div>';
          } else {
            html += '<p><strong>' + escapeHTML(s) + '</strong></p>';
          }
        });
        html += '</div></details>';
      }

      // LaTeX
      var latex = generateLaTeX(m, xs, ys, xVal, res);
      html += '<details class="latex-details"><summary>LaTeX Export</summary>';
      html += '<div class="latex-content">';
      html += '<pre class="latex-code">' + escapeHTML(latex) + '</pre>';
      html += '<button class="btn btn-small btn-copy" onclick="navigator.clipboard.writeText(this.previousElementSibling.textContent).then(function(){})">Copy LaTeX</button>';
      html += '</div></details>';

      // Share
      var shareURL = buildShareURL(xs, ys, xVal, m);
      html += '<div class="share-section">';
      html += '<input type="text" class="share-url" value="' + escapeHTML(shareURL) + '" readonly>';
      html += '<button class="btn btn-small btn-share" onclick="navigator.clipboard.writeText(this.previousElementSibling.value)">Copy Link</button>';
      html += '<button class="btn btn-small btn-secondary" onclick="window.print()">Print</button>';
      html += '</div>';

      html += '</div>';
      return html;
    }

    var result1 = runMethod(method);
    var outputHTML = '';

    if (compareMode) {
      var method2 = document.getElementById('compare-method-select').value;
      if (method2 === method) {
        // Pick a different default
        var methods = ['lagrange', 'newton-divided', 'newton-forward', 'newton-backward', 'stirling'];
        method2 = methods.find(function(m) { return m !== method && (spacing.equal || needsEqual.indexOf(m) === -1); }) || 'lagrange';
      }
      var result2 = runMethod(method2);
      outputHTML += '<div class="compare-wrapper">';
      outputHTML += buildResultHTML(method, result1);
      outputHTML += buildResultHTML(method2, result2);
      outputHTML += '</div>';
    } else {
      outputHTML = buildResultHTML(method, result1);
    }

    output.innerHTML = outputHTML;

    // Plot
    if (!result1.error) {
      var plotSection = document.getElementById('plot-section');
      plotSection.style.display = 'block';
      var canvas = document.getElementById('plot-canvas');

      // Build polyEval function using Lagrange
      var polyCoeffs = lagrangePolynomial(xs, ys);
      function polyEval(x) {
        var val = 0;
        for (var i = 0; i < polyCoeffs.length; i++) {
          val += polyCoeffs[i] * Math.pow(x, i);
        }
        return val;
      }

      plotInterpolation(canvas, xs, ys, xVal, result1.value, polyEval);
    }

    // Save to history
    if (!result1.error) {
      saveHistory({
        method: method,
        xs: xs,
        ys: ys,
        xVal: xVal,
        result: result1.value,
        n: xs.length,
        date: new Date().toLocaleString()
      });
    }
  }

  // ===== EXAMPLES =====
  function loadExample() {
    var method = document.getElementById('method-select').value;
    var exXs, exYs, exX;

    if (method === 'stirling' || (method === 'auto')) {
      // 5 equally spaced points
      document.getElementById('num-points').value = '5';
      buildDataTable(5);
      exXs = [1, 2, 3, 4, 5];
      exYs = [1, 8, 27, 64, 125];
      exX = 2.5;
    } else if (method === 'newton-forward' || method === 'newton-backward') {
      document.getElementById('num-points').value = '5';
      buildDataTable(5);
      exXs = [0, 1, 2, 3, 4];
      exYs = [1, 1, 2, 6, 24];
      exX = method === 'newton-forward' ? 0.5 : 3.5;
    } else if (method === 'lagrange' || method === 'newton-divided') {
      document.getElementById('num-points').value = '4';
      buildDataTable(4);
      exXs = [1, 2, 4, 7];
      exYs = [1, 4, 16, 49];
      exX = 3;
    } else {
      document.getElementById('num-points').value = '4';
      buildDataTable(4);
      exXs = [0, 1, 3, 4];
      exYs = [1, 3, 13, 25];
      exX = 2;
    }

    var xInputs = document.querySelectorAll('.x-input');
    var yInputs = document.querySelectorAll('.y-input');
    for (var i = 0; i < exXs.length; i++) {
      if (xInputs[i]) xInputs[i].value = exXs[i];
      if (yInputs[i]) yInputs[i].value = exYs[i];
    }
    document.getElementById('interp-x').value = exX;
  }

  // ===== IMPORT =====
  function importData() {
    var text = document.getElementById('import-text').value.trim();
    if (!text) return;

    var lines = text.split(/\n/).filter(function(l) { return l.trim(); });
    var xs = [], ys = [];
    for (var i = 0; i < lines.length; i++) {
      var parts = lines[i].trim().split(/[\s,\t]+/);
      if (parts.length >= 2) {
        xs.push(parseFloat(parts[0]));
        ys.push(parseFloat(parts[1]));
      }
    }

    if (xs.length < 2) {
      alert('Need at least 2 data points. Format: one "x y" pair per line.');
      return;
    }

    var n = Math.min(xs.length, 8);
    document.getElementById('num-points').value = n;
    buildDataTable(n);

    var xInputs = document.querySelectorAll('.x-input');
    var yInputs = document.querySelectorAll('.y-input');
    for (var i = 0; i < n; i++) {
      if (xInputs[i]) xInputs[i].value = xs[i];
      if (yInputs[i]) yInputs[i].value = ys[i];
    }

    document.getElementById('import-section').style.display = 'none';
  }

  // ===== INIT =====
  function init() {
    // Theme
    var savedTheme = localStorage.getItem('nc-interp-theme') || 'light';
    document.documentElement.setAttribute('data-theme', savedTheme);
    updateThemeBtn();

    document.getElementById('theme-toggle').addEventListener('click', function() {
      var cur = document.documentElement.getAttribute('data-theme');
      var next = cur === 'dark' ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', next);
      localStorage.setItem('nc-interp-theme', next);
      updateThemeBtn();
    });

    function updateThemeBtn() {
      var btn = document.getElementById('theme-toggle');
      btn.textContent = document.documentElement.getAttribute('data-theme') === 'dark' ? '☀️' : '🌙';
    }

    // Num points change
    document.getElementById('num-points').addEventListener('change', function() {
      buildDataTable(parseInt(this.value));
    });

    // Compare toggle
    document.getElementById('compare-toggle').addEventListener('change', function() {
      document.getElementById('compare-method-wrapper').style.display = this.checked ? 'flex' : 'none';
    });

    // Method change — show notes
    document.getElementById('method-select').addEventListener('change', function() {
      var note = document.getElementById('spacing-note');
      var m = this.value;
      if (['newton-forward', 'newton-backward', 'stirling'].indexOf(m) !== -1) {
        note.style.display = 'block';
        note.textContent = methodName(m) + ' requires equally spaced data points.';
        if (m === 'stirling') {
          note.textContent += ' Also requires an odd number of points.';
        }
      } else {
        note.style.display = 'none';
      }
    });

    // Buttons
    document.getElementById('calculate-btn').addEventListener('click', calculate);
    document.getElementById('example-btn').addEventListener('click', function() {
      loadExample();
      calculate();
    });
    document.getElementById('clear-btn').addEventListener('click', function() {
      document.getElementById('output').innerHTML = '';
      document.getElementById('plot-section').style.display = 'none';
      document.querySelectorAll('.data-input').forEach(function(inp) { inp.value = ''; inp.classList.remove('input-error'); });
      document.getElementById('interp-x').value = '';
      document.getElementById('interp-x').classList.remove('input-error');
      document.getElementById('spacing-note').style.display = 'none';
    });

    // Import
    document.getElementById('import-toggle-btn').addEventListener('click', function() {
      var sec = document.getElementById('import-section');
      sec.style.display = sec.style.display === 'none' ? 'block' : 'none';
    });
    document.getElementById('import-btn').addEventListener('click', importData);

    // History
    document.getElementById('clear-history-btn').addEventListener('click', clearHistory);
    document.getElementById('history-toggle-btn').addEventListener('click', function() {
      document.getElementById('history-panel').classList.toggle('open');
    });

    // Keyboard shortcut
    document.addEventListener('keydown', function(e) {
      if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
        e.preventDefault();
        calculate();
      }
    });

    // Build initial table
    buildDataTable(4);

    // Try loading from URL
    var loaded = loadFromURL();
    if (loaded) {
      setTimeout(calculate, 100);
    }

    // Render history
    renderHistory();
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }

})();
