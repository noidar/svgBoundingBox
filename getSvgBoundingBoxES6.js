class BoundingBox {
    constructor(x1, y1, x2, y2) {
        this.x1 = Number.NaN;
        this.y1 = Number.NaN;
        this.x2 = Number.NaN;
        this.y2 = Number.NaN;
      
        this.addPoint(x1, y1);
        this.addPoint(x2, y2);
      }
  
    width() {
      return this.x2 - this.x1;
    }
  
    height() {
      return this.y2 - this.y1;
    }
  
    addPoint(x, y) {
      if (x != null) {
        if (isNaN(this.x1) || isNaN(this.x2)) {
          this.x1 = x;
          this.x2 = x;
        }
        if (x < this.x1) this.x1 = x;
        if (x > this.x2) this.x2 = x;
      }
  
      if (y != null) {
        if (isNaN(this.y1) || isNaN(this.y2)) {
          this.y1 = y;
          this.y2 = y;
        }
        if (y < this.y1) this.y1 = y;
        if (y > this.y2) this.y2 = y;
      }
    }
  
    addX(x) {
      this.addPoint(x, null);
    }
  
    addY(y) {
      this.addPoint(null, y);
    }
  
    addQuadraticCurve(p0x, p0y, p1x, p1y, p2x, p2y) {
      const cp1x = p0x + 2 / 3 * (p1x - p0x); // CP1 = QP0 + 2/3 *(QP1-QP0)
      const cp1y = p0y + 2 / 3 * (p1y - p0y); // CP1 = QP0 + 2/3 *(QP1-QP0)
      const cp2x = cp1x + 1 / 3 * (p2x - p0x); // CP2 = CP1 + 1/3 *(QP2-QP0)
      const cp2y = cp1y + 1 / 3 * (p2y - p0y); // CP2 = CP1 + 1/3 *(QP2-QP0)
      this.addBezierCurve(p0x, p0y, cp1x, cp2x, cp1y, cp2y, p2x, p2y);
    }
  
    addBezierCurve(p0x, p0y, p1x, p1y, p2x, p2y, p3x, p3y) {
      // from http://blog.hackers-cafe.net/2009/06/how-to-calculate-bezier-curves-bounding.html
      let i;
  
      const p0 = [p0x, p0y];
      const p1 = [p1x, p1y];
      const p2 = [p2x, p2y];
      const p3 = [p3x, p3y];
  
      this.addPoint(p0[0], p0[1]);
      this.addPoint(p3[0], p3[1]);
  
      for (i = 0; i <= 1; i++) {
        const f = t => (1 - t) ** 3 * p0[i]
          + 3 * ((1 - t) ** 2) * t * p1[i]
          + 3 * (1 - t) * (t ** 2) * p2[i]
          + t ** 3 * p3[i];
  
        const b = 6 * p0[i] - 12 * p1[i] + 6 * p2[i];
        const a = -3 * p0[i] + 9 * p1[i] - 9 * p2[i] + 3 * p3[i];
        const c = 3 * p1[i] - 3 * p0[i];
  
        if (a == 0) {
          if (b == 0) continue;
          const t = -c / b;
          if (0 < t && t < 1) {
            if (i == 0) this.addX(f(t));
            if (i == 1) this.addY(f(t));
          }
          continue;
        }
  
        const b2ac = b ** 2 - 4 * c * a;
        if (b2ac < 0) continue;
        const t1 = (-b + Math.sqrt(b2ac)) / (2 * a);
        if (0 < t1 && t1 < 1) {
          if (i == 0) this.addX(f(t1));
          if (i == 1) this.addY(f(t1));
        }
        const t2 = (-b - Math.sqrt(b2ac)) / (2 * a);
        if (0 < t2 && t2 < 1) {
          if (i == 0) this.addX(f(t2));
          if (i == 1) this.addY(f(t2));
        }
      }
    }
  }
  
  
  
  function BoundingBoxView(boundingBox) {
  
    this.x1 = this.minX = boundingBox.x1 || 0;
    this.y1 = this.minY = boundingBox.y1 || 0;
    this.x2 = this.maxX = boundingBox.x2 || 0;
    this.y2 = this.maxY = boundingBox.y2 || 0;
    this.width = boundingBox.width() || 0;
    this.height = boundingBox.height() || 0;
  
  }
  
  const paramCounts = { a: 7, c: 6, h: 1, l: 2, m: 2, r: 4, q: 4, s: 4, t: 2, v: 1, z: 0 };
  
  function isSpace(ch) {
  
    return /^[ \t\r\n]$/.test(ch);
  }
  
  function isCommand(ch) {
  
    return /^[mzlhvcsqtar]$/i.test(ch);
  }
  
  function isDigit(ch) {
  
    return /^\d$/.test(ch);
  }
  
  function isDigitStart(ch) {
  
    return /^[\d+\-\.]$/.test(ch);
  }
  
  function isSign(ch) {
  
    return /^[+-]$/.test(ch);
  }
  
  function skipSpaces(state) {
    while (state.index < state.max && isSpace(state.path.charAt(state.index))) {
      state.index++;
    }
  }
  
  function scanParam(state) {
    const start = state.index;
    let index = start;
    const max = state.max;
    let zeroFirst = false;
    let hasCeiling = false;
    let hasDecimal = false;
    let hasDot = false;
    let ch;
  
    if (index >= max) {
      state.err = `SvgPath: missed param (at pos ${index})`;
      return;
    }
    ch = state.path.charCodeAt(index);
  
    if (isSign(state.path.charAt(index))) {
      index++;
      ch = (index < max) ? state.path.charCodeAt(index) : 0;
    }
  
    // This logic is shamelessly borrowed from Esprima
    // https://github.com/ariya/esprimas
    //
    if (!isDigit(String.fromCharCode(ch)) && ch !== 0x2E/* . */) {
      state.err = `SvgPath: param should start with 0..9 or \`.\` (at pos ${index})`;
      return;
    }
  
    if (ch !== 0x2E/* . */) {
      zeroFirst = (ch === 0x30/* 0 */);
      index++;
  
      ch = (index < max) ? state.path.charCodeAt(index) : 0;
  
      if (zeroFirst && index < max) {
        // decimal number starts with '0' such as '09' is illegal.
        if (ch && isDigit(String.fromCharCode(ch))) {
          state.err = `SvgPath: numbers started with \`0\` such as \`09\` are ilegal (at pos ${start})`;
          return;
        }
      }
  
      while (index < max && isDigit(state.path.charAt(index))) {
        index++;
        hasCeiling = true;
      }
      ch = (index < max) ? state.path.charCodeAt(index) : 0;
    }
  
    if (ch === 0x2E/* . */) {
      hasDot = true;
      index++;
      while (isDigit(state.path.charAt(index))) {
        index++;
        hasDecimal = true;
      }
      ch = (index < max) ? state.path.charCodeAt(index) : 0;
    }
  
    if (ch === 0x65/* e */ || ch === 0x45/* E */) {
      if (hasDot && !hasCeiling && !hasDecimal) {
        state.err = `SvgPath: invalid float exponent (at pos ${index})`;
        return;
      }
  
      index++;
  
      ch = (index < max) ? state.path.charCodeAt(index) : 0;
      if (ch === 0x2B/* + */ || ch === 0x2D/* - */) {
        index++;
      }
      if (index < max && isDigit(state.path.charAt(index))) {
        while (index < max && isDigit(state.path.charAt(index))) {
          index++;
        }
      } else {
        state.err = `SvgPath: invalid float exponent (at pos ${index})`;
        return;
      }
    }
  
    state.index = index;
    state.param = parseFloat(state.path.slice(start, index)) + 0.0;
  }
  
  function finalizeSegment(state) {
    let cmd;
    let cmdLC;
  
    // Process duplicated commands (without comand name)
  
    // This logic is shamelessly borrowed from Raphael
    // https://github.com/DmitryBaranovskiy/raphael/
    //
    cmd   = state.path[state.segmentStart];
    cmdLC = cmd.toLowerCase();
  
    let params = state.data;
  
    if (cmdLC === 'm' && params.length > 2) {
      state.result.push([ cmd, params[0], params[1] ]);
      params = params.slice(2);
      cmdLC = 'l';
      cmd = (cmd === 'm') ? 'l' : 'L';
    }
  
    if (cmdLC === 'r') {
      state.result.push([ cmd ].concat(params));
    } else {
  
      while (params.length >= paramCounts[cmdLC]) {
        state.result.push([ cmd ].concat(params.splice(0, paramCounts[cmdLC])));
        if (!paramCounts[cmdLC]) {
          break;
        }
      }
    }
  }
  
  
  function scanSegment(state) {
    const max = state.max;
    let comma_found;
    let need_params;
    let i;
  
    state.segmentStart = state.index;
  
    if (!isCommand(state.path.charAt(state.index))) {
      state.err = `SvgPath: bad command ${state.path[state.index]} (at pos ${state.index})`;
      return;
    }
  
    need_params = paramCounts[state.path[state.index].toLowerCase()];
  
    state.index++;
    skipSpaces(state);
  
    state.data = [];
  
    if (!need_params) {
      // Z
      finalizeSegment(state);
      return;
    }
  
    comma_found = false;
  
    for (;;) {
      for (i = need_params; i > 0; i--) {
        scanParam(state);
        if (state.err.length) {
          return;
        }
        state.data.push(state.param);
  
        skipSpaces(state);
        comma_found = false;
  
        if (state.index < max && state.path.charCodeAt(state.index) === 0x2C/* , */) {
          state.index++;
          skipSpaces(state);
          comma_found = true;
        }
      }
  
      // after ',' param is mandatory
      if (comma_found) {
        continue;
      }
  
      if (state.index >= state.max) {
        break;
      }
  
      // Stop on next segment
      if (!isDigitStart(state.path.charAt(state.index))) {
        break;
      }
    }
  
    finalizeSegment(state);
  }
  
  
  /* Returns array of segments:
   *
   * [
   *   [ command, coord1, coord2, ... ]
   * ]
   */
  function pathParse(svgPath) {
    const state = {
      index: 0,
      path: svgPath,
      max: svgPath.length,
      result: [],
      param: 0.0,
      err: '',
      segmentStart: 0,
      data: []
    };
    const max = state.max;
  
    skipSpaces(state);
  
    while (state.index < max && !state.err.length) {
      scanSegment(state);
    }
  
    if (state.err.length) {
      state.result = [];
  
    } else if (state.result.length) {
  
      if (!'mM'.includes(state.result[0][0])) {
        state.err = 'SvgPath: string should start with `M` or `m`';
        state.result = [];
      } else {
        state.result[0][0] = 'M';
      }
    }
  
    return {
      err: state.err,
      segments: state.result
    };
  }
  
  const TAU = Math.PI * 2;
  
  
  /* eslint-disable space-infix-ops */
  
  // Calculate an angle between two vectors
  //
  function vector_angle(ux, uy, vx, vy) {
    const sign = (ux * vy - uy * vx < 0) ? -1 : 1;
    const umag = Math.sqrt(ux * ux + uy * uy);
    const vmag = Math.sqrt(ux * ux + uy * uy);
    const dot  = ux * vx + uy * vy;
    let div  = dot / (umag * vmag);
  
    // rounding errors, e.g. -1.0000000000000002 can screw up this
    if (div > 1) { div = 1; }
    if (div < -1) { div = -1; }
  
    return sign * Math.acos(div);
  }
  
  
  // Convert from endpoint to center parameterization,
  // see http://www.w3.org/TR/SVG11/implnote.html#ArcImplementationNotes
  //
  // Return [cx, cy, theta1, delta_theta]
  //
  function get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi) {
    // Step 1.
    //
    // Moving an ellipse so origin will be the middlepoint between our two
    // points. After that, rotate it to line up ellipse axes with coordinate
    // axes.
    //
    const x1p =  cos_phi*(x1-x2)/2 + sin_phi*(y1-y2)/2;
    const y1p = -sin_phi*(x1-x2)/2 + cos_phi*(y1-y2)/2;
  
    const rx_sq  =  rx * rx;
    const ry_sq  =  ry * ry;
    const x1p_sq = x1p * x1p;
    const y1p_sq = y1p * y1p;
  
    // Step 2.
    //
    // Compute coordinates of the centre of this ellipse (cx', cy')
    // in the new coordinate system.
    //
    let radicant = (rx_sq * ry_sq) - (rx_sq * y1p_sq) - (ry_sq * x1p_sq);
  
    if (radicant < 0) {
      // due to rounding errors it might be e.g. -1.3877787807814457e-17
      radicant = 0;
    }
  
    radicant /=   (rx_sq * y1p_sq) + (ry_sq * x1p_sq);
    radicant = Math.sqrt(radicant) * (fa === fs ? -1 : 1);
  
    const cxp = radicant *  rx/ry * y1p;
    const cyp = radicant * -ry/rx * x1p;
  
    // Step 3.
    //
    // Transform back to get centre coordinates (cx, cy) in the original
    // coordinate system.
    //
    const cx = cos_phi*cxp - sin_phi*cyp + (x1+x2)/2;
    const cy = sin_phi*cxp + cos_phi*cyp + (y1+y2)/2;
  
    // Step 4.
    //
    // Compute angles (theta1, delta_theta).
    //
    const v1x =  (x1p - cxp) / rx;
    const v1y =  (y1p - cyp) / ry;
    const v2x = (-x1p - cxp) / rx;
    const v2y = (-y1p - cyp) / ry;
  
    const theta1 = vector_angle(1, 0, v1x, v1y);
    let delta_theta = vector_angle(v1x, v1y, v2x, v2y);
  
    if (fs === 0 && delta_theta > 0) {
      delta_theta -= TAU;
    }
    if (fs === 1 && delta_theta < 0) {
      delta_theta += TAU;
    }
  
    return [ cx, cy, theta1, delta_theta ];
  }
  
  //
  // Approximate one unit arc segment with bézier curves,
  // see http://math.stackexchange.com/questions/873224
  //
  function approximate_unit_arc(theta1, delta_theta) {
    const alpha = 4/3 * Math.tan(delta_theta/4);
  
    const x1 = Math.cos(theta1);
    const y1 = Math.sin(theta1);
    const x2 = Math.cos(theta1 + delta_theta);
    const y2 = Math.sin(theta1 + delta_theta);
  
    return [ x1, y1, x1 - y1*alpha, y1 + x1*alpha, x2 + y2*alpha, y2 - x2*alpha, x2, y2 ];
  }
  
  function a2c(x1, y1, x2, y2, fa, fs, rx, ry, phi) {
    const sin_phi = Math.sin(phi * TAU / 360);
    const cos_phi = Math.cos(phi * TAU / 360);
  
    // Make sure radii are valid
    //
    const x1p =  cos_phi*(x1-x2)/2 + sin_phi*(y1-y2)/2;
    const y1p = -sin_phi*(x1-x2)/2 + cos_phi*(y1-y2)/2;
  
    if (x1p === 0 && y1p === 0) {
      // we're asked to draw line to itself
      return [];
    }
  
    if (rx === 0 || ry === 0) {
      // one of the radii is zero
      return [];
    }
  
  
    // Compensate out-of-range radii
    //
    rx = Math.abs(rx);
    ry = Math.abs(ry);
  
    const lambda = (x1p * x1p) / (rx * rx) + (y1p * y1p) / (ry * ry);
    if (lambda > 1) {
      rx *= Math.sqrt(lambda);
      ry *= Math.sqrt(lambda);
    }
  
  
    // Get center parameters (cx, cy, theta1, delta_theta)
    //
    const cc = get_arc_center(x1, y1, x2, y2, fa, fs, rx, ry, sin_phi, cos_phi);
  
    const result = [];
    let theta1 = cc[2];
    let delta_theta = cc[3];
  
    // Split an arc to multiple segments, so each segment
    // will be less than τ/4 (= 90°)
    //
    const segments = Math.max(Math.ceil(Math.abs(delta_theta) / (TAU / 4)), 1);
    delta_theta /= segments;
  
    for (let i = 0; i < segments; i++) {
      result.push(approximate_unit_arc(theta1, delta_theta));
      theta1 += delta_theta;
    }
  
    // We have a bezier approximation of a unit circle,
    // now need to transform back to the original ellipse
    //
    return result.map(curve => {
      for (let i = 0; i < curve.length; i += 2) {
        let x = curve[i + 0];
        let y = curve[i + 1];
  
        // scale
        x *= rx;
        y *= ry;
  
        // rotate
        const xp = cos_phi*x - sin_phi*y;
        const yp = sin_phi*x + cos_phi*y;
  
        // translate
        curve[i + 0] = xp + cc[0];
        curve[i + 1] = yp + cc[1];
      }
  
      return curve;
    });
  }
  
  //#################################################
  // Class constructor
  //*********************************************************************** */
  class SvgPath {
    constructor(path) {
    
      const pstate = pathParse(path);
    
      // Array of path segments.
      // Each segment is array [command, param1, param2, ...]
      this.segments = pstate.segments;
    
      // Error message on parse error.
      this.err      = pstate.err;
    
      // Transforms stack for lazy evaluation
      this.__stack    = [];
    }
  
    // Apply stacked commands
    //
    __evaluateStack() {
      let m;
      let i;
  
      if (!this.__stack.length) { return; }
  
      if (this.__stack.length === 1) {
        this.__matrix(this.__stack[0]);
        this.__stack = [];
        return;
      }
  
      m = matrix();
      i = this.__stack.length;
  
      while (--i >= 0) {
        m.matrix(this.__stack[i].toArray());
      }
  
      this.__matrix(m);
      this.__stack = [];
    }
  
    // Apply iterator function to all segments. If function returns result,
    // current segment will be replaced to array of returned segments.
    // If empty array is returned, current regment will be deleted.
    //
    iterate(iterator, keepLazyStack) {
      const segments = this.segments;
      const replacements = {};
      let needReplace = false;
      let lastX = 0;
      let lastY = 0;
      let countourStartX = 0;
      let countourStartY = 0;
      let i;
      let j;
      let newSegments;
  
      if (!keepLazyStack) {
        this.__evaluateStack();
      }
  
      segments.forEach((s, index) => {
    
        const res = iterator(s, index, lastX, lastY);
    
        if (Array.isArray(res)) {
          replacements[index] = res;
          needReplace = true;
        }
    
        const isRelative = (s[0] === s[0].toLowerCase());
    
        // calculate absolute X and Y
        switch (s[0]) {
          case 'm':
          case 'M':
            lastX = s[1] + (isRelative ? lastX : 0);
            lastY = s[2] + (isRelative ? lastY : 0);
            countourStartX = lastX;
            countourStartY = lastY;
            return;
    
          case 'h':
          case 'H':
            lastX = s[1] + (isRelative ? lastX : 0);
            return;
    
          case 'v':
          case 'V':
            lastY = s[1] + (isRelative ? lastY : 0);
            return;
    
          case 'z':
          case 'Z':
            // That make sence for multiple contours
            lastX = countourStartX;
            lastY = countourStartY;
            return;
    
          default:
            lastX = s[s.length - 2] + (isRelative ? lastX : 0);
            lastY = s[s.length - 1] + (isRelative ? lastY : 0);
        }
      });
  
      // Replace segments if iterator return results
  
      if (!needReplace) { return this; }
  
      newSegments = [];
  
      for (i = 0; i < segments.length; i++) {
        if (typeof replacements[i] !== 'undefined') {
          for (j = 0; j < replacements[i].length; j++) {
            newSegments.push(replacements[i][j]);
          }
        } else {
          newSegments.push(segments[i]);
        }
      }
  
      this.segments = newSegments;
  
      return this;
    }
  
    // Converts segments from relative to absolute
    //
    abs() {
    
      this.iterate((s, index, x, y) => {
        const name = s[0];
        const nameUC = name.toUpperCase();
        let i;
  
        // Skip absolute commands
        if (name === nameUC) { return; }
  
        s[0] = nameUC;
  
        switch (name) {
          case 'v':
            // v has shifted coords parity
            s[1] += y;
            return;
    
          case 'a':
            // ARC is: ['A', rx, ry, x-axis-rotation, large-arc-flag, sweep-flag, x, y]
            // touch x, y only
            s[6] += x;
            s[7] += y;
            return;
    
          default:
            for (i = 1; i < s.length; i++) {
              s[i] += i % 2 ? x : y; // odd values are X, even - Y
            }
        }
      }, true);
    
      return this;
    }
  
    // Converts arcs to cubic bézier curves
    //
    unarc() {
      this.iterate((s, index, x, y) => {
        let new_segments;
        let nextX;
        let nextY;
        const result = [];
        const name = s[0];
  
        // Skip anything except arcs
        if (name !== 'A' && name !== 'a') { return null; }
  
        if (name === 'a') {
          // convert relative arc coordinates to absolute
          nextX = x + s[6];
          nextY = y + s[7];
        } else {
          nextX = s[6];
          nextY = s[7];
        }
  
        new_segments = a2c(x, y, nextX, nextY, s[4], s[5], s[1], s[2], s[3]);
  
        // Degenerated arcs can be ignored by renderer, but should not be dropped
        // to avoid collisions with `S A S` and so on. Replace with empty line.
        if (new_segments.length === 0) {
          return [ [ s[0] === 'a' ? 'l' : 'L', s[6], s[7] ] ];
        }
  
        new_segments.forEach(s => {
          result.push([ 'C', s[2], s[3], s[4], s[5], s[6], s[7] ]);
        });
  
        return result;
      });
    
      return this;
    }
  
    // Converts smooth curves (with missed control point) to generic curves
    //
    unshort() {
      const segments = this.segments;
      let prevControlX;
      let prevControlY;
      let prevSegment;
      let curControlX;
      let curControlY;
  
      // TODO: add lazy evaluation flag when relative commands supported
  
      this.iterate((s, idx, x, y) => {
        const name = s[0];
        const nameUC = name.toUpperCase();
        let isRelative;
  
        // First command MUST be M|m, it's safe to skip.
        // Protect from access to [-1] for sure.
        if (!idx) { return; }
  
        if (nameUC === 'T') { // quadratic curve
          isRelative = (name === 't');
    
          prevSegment = segments[idx - 1];
    
          if (prevSegment[0] === 'Q') {
            prevControlX = prevSegment[1] - x;
            prevControlY = prevSegment[2] - y;
          } else if (prevSegment[0] === 'q') {
            prevControlX = prevSegment[1] - prevSegment[3];
            prevControlY = prevSegment[2] - prevSegment[4];
          } else {
            prevControlX = 0;
            prevControlY = 0;
          }
    
          curControlX = -prevControlX;
          curControlY = -prevControlY;
    
          if (!isRelative) {
            curControlX += x;
            curControlY += y;
          }
    
          segments[idx] = [
            isRelative ? 'q' : 'Q',
            curControlX, curControlY,
            s[1], s[2]
          ];
    
        } else if (nameUC === 'S') { // cubic curve
          isRelative = (name === 's');
    
          prevSegment = segments[idx - 1];
    
          if (prevSegment[0] === 'C') {
            prevControlX = prevSegment[3] - x;
            prevControlY = prevSegment[4] - y;
          } else if (prevSegment[0] === 'c') {
            prevControlX = prevSegment[3] - prevSegment[5];
            prevControlY = prevSegment[4] - prevSegment[6];
          } else {
            prevControlX = 0;
            prevControlY = 0;
          }
    
          curControlX = -prevControlX;
          curControlY = -prevControlY;
    
          if (!isRelative) {
            curControlX += x;
            curControlY += y;
          }
    
          segments[idx] = [
            isRelative ? 'c' : 'C',
            curControlX, curControlY,
            s[1], s[2], s[3], s[4]
          ];
        }
      });
  
      return this;
    }
  }
  

  class Path {
    constructor(d) {
      this.d = d;
    }
  
    getBoundingBox() {
      let pathDriver;
      let boundingBox;
  
    pathDriver = new SvgPath(this.d);
      boundingBox = new BoundingBox();
      pathDriver
        .abs()
        .unarc()
        .unshort()
        .iterate((seg, index, x, y) => {
  
          switch(seg[0]) {
            case 'M':
            case 'L':
              boundingBox.addPoint(
                seg[1],
                seg[2]
              );
              break;
            case 'H':
              boundingBox.addX(seg[1]);
              break;
            case 'V':
              boundingBox.addY(seg[1]);
              break;
            case 'Q':
              boundingBox.addQuadraticCurve(
                x,
                y,
                seg[1],
                seg[2],
                seg[3],
                seg[4]
              );
              break;
            case 'C':
              boundingBox.addBezierCurve(
                x,
                y,
                seg[1],
                seg[2],
                seg[3],
                seg[4],
                seg[5],
                seg[6]
              );
              break;
          }
  
        });
  
      return new BoundingBoxView(boundingBox);
    }
  }
  
const response = new Path(path).getBoundingBox();

console.table(response)
