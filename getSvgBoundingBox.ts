class BoundingBox {

	public x1: number;
	public y1: number;
	public x2: number;
	public y2: number;

	constructor(x1 = Number.NaN, y1 = Number.NaN, x2 = Number.NaN, y2 = Number.NaN) {
		this.x1 = x1;
		this.y1 = y1;
		this.x2 = x2;
		this.y2 = y2;

		this.addPoint(x1, y1);
		this.addPoint(x2, y2);
	}

	public width() {
		return this.x2 - this.x1;
	}

	public height() {
		return this.y2 - this.y1;
	}

	public addPoint(x: number, y: number) {
		if (x !== null) {
			if (isNaN(this.x1) || isNaN(this.x2)) {
				this.x1 = x;
				this.x2 = x;
			}
			if (x < this.x1) { this.x1 = x; }
			if (x > this.x2) { this.x2 = x; }
		}

		if (y !== null) {
			if (isNaN(this.y1) || isNaN(this.y2)) {
				this.y1 = y;
				this.y2 = y;
			}
			if (y < this.y1) { this.y1 = y; }
			if (y > this.y2) { this.y2 = y; }
		}
	}

	public addX(x: number) {
		this.addPoint(x, null);
	}

	public addY(y: number) {
		this.addPoint(null, y);
	}

	public addQuadraticCurve(p0x: number, p0y: number, p1x: number, p1y: number, p2x: number, p2y: number) {
		const cp1x = p0x + 2 / 3 * (p1x - p0x); // CP1 = QP0 + 2/3 *(QP1-QP0)
		const cp1y = p0y + 2 / 3 * (p1y - p0y); // CP1 = QP0 + 2/3 *(QP1-QP0)
		const cp2x = cp1x + 1 / 3 * (p2x - p0x); // CP2 = CP1 + 1/3 *(QP2-QP0)
		const cp2y = cp1y + 1 / 3 * (p2y - p0y); // CP2 = CP1 + 1/3 *(QP2-QP0)
		this.addBezierCurve(p0x, p0y, cp1x, cp2x, cp1y, cp2y, p2x, p2y);
	}

	public addBezierCurve(p0x: number, p0y: number, p1x: number, p1y: number, p2x: number, p2y: number, p3x: number, p3y: number) {
		// from http://blog.hackers-cafe.net/2009/06/how-to-calculate-bezier-curves-bounding.html

		const p0 = [p0x, p0y];
		const p1 = [p1x, p1y];
		const p2 = [p2x, p2y];
		const p3 = [p3x, p3y];

		this.addPoint(p0[0], p0[1]);
		this.addPoint(p3[0], p3[1]);

		for (let i = 0; i <= 1; i++) {
			const f = (t) => (1 - t) ** 3 * p0[i]
				+ 3 * ((1 - t) ** 2) * t * p1[i]
				+ 3 * (1 - t) * (t ** 2) * p2[i]
				+ t ** 3 * p3[i];

			const b = 6 * p0[i] - 12 * p1[i] + 6 * p2[i];
			const a = -3 * p0[i] + 9 * p1[i] - 9 * p2[i] + 3 * p3[i];
			const c = 3 * p1[i] - 3 * p0[i];

			if (a === 0) {
				if (b === 0) { continue; }
				const t = -c / b;
				if (0 < t && t < 1) {
					if (i === 0) { this.addX(f(t)); }
					if (i === 1) { this.addY(f(t)); }
				}
				continue;
			}

			const b2ac = b ** 2 - 4 * c * a;
			if (b2ac < 0) { continue; }
			const t1 = (-b + Math.sqrt(b2ac)) / (2 * a);
			if (0 < t1 && t1 < 1) {
				if (i === 0) { this.addX(f(t1)); }
				if (i === 1) { this.addY(f(t1)); }
			}
			const t2 = (-b - Math.sqrt(b2ac)) / (2 * a);
			if (0 < t2 && t2 < 1) {
				if (i === 0) { this.addX(f(t2)); }
				if (i === 1) { this.addY(f(t2)); }
			}
		}
	}
}

class SvgPath {
	private segments: any [][];
	private err: string;
	private stack;

	constructor(path: string) {

		const pstate = Path.pathParse(path);

		// Array of path segments.
		// Each segment is array [command, param1, param2, ...]
		this.segments = pstate.segments;

		// Error message on parse error.
		this.err = pstate.err;

		// Transforms stack for lazy evaluation
		this.stack = [];
	}

	// Apply stacked commands

	public evaluateStack() {

		if (!this.stack.length) { return; }

		if (this.stack.length === 1) {
			this.stack = [];
			return;
		}

		this.stack = [];
	}

	// Apply iterator function to all segments. If function returns result,
	// current segment will be replaced to array of returned segments.
	// If empty array is returned, current regment will be deleted.
	//
	public iterate(iterator, keepLazyStack = null) {
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

		if (!keepLazyStack && keepLazyStack !== null) {
			this.evaluateStack();
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
	public abs() {

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
	public unarc() {
		this.iterate((s, index, x, y) => {
			let newSegments;
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

			newSegments = Path.a2c(x, y, nextX, nextY, s[4], s[5], s[1], s[2], s[3]);

			// Degenerated arcs can be ignored by renderer, but should not be dropped
			// to avoid collisions with `S A S` and so on. Replace with empty line.
			if (newSegments.length === 0) {
				return [[s[0] === 'a' ? 'l' : 'L', s[6], s[7]]];
			}

			newSegments.forEach((i) => {
				result.push(['C', i[2], i[3], i[4], i[5], i[6], i[7]]);
			});

			return result;
		});

		return this;
	}

	// Converts smooth curves (with missed control point) to generic curves
	//
	public unshort() {
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

export class Path {
	private d;
	private static paramCounts = { a: 7, c: 6, h: 1, l: 2, m: 2, r: 4, q: 4, s: 4, t: 2, v: 1, z: 0 };
	private static readonly TAU = Math.PI * 2;

	constructor(d) {
		this.d = d;
	}


	
	public static getBoudingBox(path) {
		return new Path(path).getBoundingBox();
	}

	private getBoundingBox() {

		const pathDriver = new SvgPath(this.d);
		const boundingBox = new BoundingBox();

		pathDriver
			.abs()
			.unarc()
			.unshort()
			.iterate((seg, index, x, y) => {

				switch (seg[0]) {
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

		return {
			minX: boundingBox.x1 || 0,
			minY: boundingBox.y1 || 0,
			maxX: boundingBox.x2 || 0,
			maxY: boundingBox.y2 || 0,
			width: boundingBox.width() || 0,
			height: boundingBox.height() || 0
		};
	}

	public static isSpace = (ch) => {
		return /^[ \t\r\n]$/.test(ch);
	}

	public static isCommand(ch) {
		return /^[mzlhvcsqtar]$/i.test(ch);
	}

	public static isDigit(ch) {
		return /^\d$/.test(ch);
	}

	public static isDigitStart(ch) {
		return /^[\d+\-\.]$/.test(ch);
	}

	public static isSign(ch) {
		return /^[+-]$/.test(ch);
	}

	public static skipSpaces(state) {
		while (state.index < state.max && this.isSpace(state.path.charAt(state.index))) {
			state.index++;
		}
	}

	public static scanParam(state) {
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

		if (this.isSign(state.path.charAt(index))) {
			index++;
			ch = (index < max) ? state.path.charCodeAt(index) : 0;
		}

		// This logic is shamelessly borrowed from Esprima
		// https://github.com/ariya/esprimas
		//
		if (!this.isDigit(String.fromCharCode(ch)) && ch !== 0x2E/* . */) {
			state.err = `SvgPath: param should start with 0..9 or \`.\` (at pos ${index})`;
			return;
		}

		if (ch !== 0x2E/* . */) {
			zeroFirst = (ch === 0x30/* 0 */);
			index++;

			ch = (index < max) ? state.path.charCodeAt(index) : 0;

			if (zeroFirst && index < max) {
				// decimal number starts with '0' such as '09' is illegal.
				if (ch && this.isDigit(String.fromCharCode(ch))) {
					state.err = `SvgPath: numbers started with \`0\` such as \`09\` are ilegal (at pos ${start})`;
					return;
				}
			}

			while (index < max && this.isDigit(state.path.charAt(index))) {
				index++;
				hasCeiling = true;
			}
			ch = (index < max) ? state.path.charCodeAt(index) : 0;
		}

		if (ch === 0x2E/* . */) {
			hasDot = true;
			index++;
			while (this.isDigit(state.path.charAt(index))) {
				index++;
				hasDecimal = true;
			}
			ch = (index < max) ? state.path.charCodeAt(index) : 0;
		}

		if (ch === 0x65 /* e */ || ch === 0x45/* E */) {
			if (hasDot && !hasCeiling && !hasDecimal) {
				state.err = `SvgPath: invalid float exponent (at pos ${index})`;
				return;
			}

			index++;

			ch = (index < max) ? state.path.charCodeAt(index) : 0;
			if (ch === 0x2B /* + */ || ch === 0x2D/* - */) {
				index++;
			}
			if (index < max && this.isDigit(state.path.charAt(index))) {
				while (index < max && this.isDigit(state.path.charAt(index))) {
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

	public static finalizeSegment(state) {
		let cmd;
		let cmdLC;

		// Process duplicated commands (without comand name)

		// This logic is shamelessly borrowed from Raphael
		// https://github.com/DmitryBaranovskiy/raphael/
		//
		cmd = state.path[state.segmentStart];
		cmdLC = cmd.toLowerCase();

		let params = state.data;

		if (cmdLC === 'm' && params.length > 2) {
			state.result.push([cmd, params[0], params[1]]);
			params = params.slice(2);
			cmdLC = 'l';
			cmd = (cmd === 'm') ? 'l' : 'L';
		}

		if (cmdLC === 'r') {
			state.result.push([cmd].concat(params));
		} else {

			while (params.length >= Path.paramCounts[cmdLC]) {
				state.result.push([cmd].concat(params.splice(0, Path.paramCounts[cmdLC])));
				if (!Path.paramCounts[cmdLC]) {
					break;
				}
			}
		}
	}

	public static scanSegment(state) {
		const max = state.max;
		let commaFound;
		let needParams;
		let i;

		state.segmentStart = state.index;

		if (!this.isCommand(state.path.charAt(state.index))) {
			state.err = `SvgPath: bad command ${state.path[state.index]} (at pos ${state.index})`;
			return;
		}

		needParams = Path.paramCounts[state.path[state.index].toLowerCase()];

		state.index++;
		this.skipSpaces(state);

		state.data = [];

		if (!needParams) {
			// Z
			this.finalizeSegment(state);
			return;
		}

		commaFound = false;

		for (; ;) {
			for (i = needParams; i > 0; i--) {
				this.scanParam(state);
				if (state.err.length) {
					return;
				}
				state.data.push(state.param);

				this.skipSpaces(state);
				commaFound = false;

				if (state.index < max && state.path.charCodeAt(state.index) === 0x2C/* , */) {
					state.index++;
					this.skipSpaces(state);
					commaFound = true;
				}
			}

			// after ',' param is mandatory
			if (commaFound) {
				continue;
			}

			if (state.index >= state.max) {
				break;
			}

			// Stop on next segment
			if (!this.isDigitStart(state.path.charAt(state.index))) {
				break;
			}
		}

		this.finalizeSegment(state);
	}

	public static pathParse(svgPath) {
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

		this.skipSpaces(state);

		while (state.index < max && !state.err.length) {
			this.scanSegment(state);
		}

		if (state.err.length) {
			state.result = [];

		} else if (state.result.length) {
			state.result[0][0] = 'M';
		}

		return {
			err: state.err,
			segments: state.result
		};
	}

	public static vectorAngle(ux, uy, vx, vy) {
		const sign = (ux * vy - uy * vx < 0) ? -1 : 1;
		const umag = Math.sqrt(ux * ux + uy * uy);
		const vmag = Math.sqrt(ux * ux + uy * uy);
		const dot = ux * vx + uy * vy;
		let div = dot / (umag * vmag);

		// rounding errors, e.g. -1.0000000000000002 can screw up this
		if (div > 1) { div = 1; }
		if (div < -1) { div = -1; }

		return sign * Math.acos(div);
	}

	public static getArcCenter(x1, y1, x2, y2, fa, fs, rx, ry, sinPhi, cosPhi) {
		// Step 1.
		//
		// Moving an ellipse so origin will be the middlepoint between our two
		// points. After that, rotate it to line up ellipse axes with coordinate
		// axes.
		//
		const x1p = cosPhi * (x1 - x2) / 2 + sinPhi * (y1 - y2) / 2;
		const y1p = -sinPhi * (x1 - x2) / 2 + cosPhi * (y1 - y2) / 2;

		const rxSq = rx * rx;
		const rySq = ry * ry;
		const x1pSq = x1p * x1p;
		const y1pSq = y1p * y1p;

		// Step 2.
		//
		// Compute coordinates of the centre of this ellipse (cx', cy')
		// in the new coordinate system.
		//
		let radicant = (rxSq * rySq) - (rxSq * y1pSq) - (rySq * x1pSq);

		if (radicant < 0) {
			// due to rounding errors it might be e.g. -1.3877787807814457e-17
			radicant = 0;
		}

		radicant /= (rxSq * y1pSq) + (rySq * x1pSq);
		radicant = Math.sqrt(radicant) * (fa === fs ? -1 : 1);

		const cxp = radicant * rx / ry * y1p;
		const cyp = radicant * -ry / rx * x1p;

		// Step 3.
		//
		// Transform back to get centre coordinates (cx, cy) in the original
		// coordinate system.
		//
		const cx = cosPhi * cxp - sinPhi * cyp + (x1 + x2) / 2;
		const cy = sinPhi * cxp + cosPhi * cyp + (y1 + y2) / 2;

		// Step 4.
		//
		// Compute angles (theta1, delta_theta).
		//
		const v1x = (x1p - cxp) / rx;
		const v1y = (y1p - cyp) / ry;
		const v2x = (-x1p - cxp) / rx;
		const v2y = (-y1p - cyp) / ry;

		const theta1 = this.vectorAngle(1, 0, v1x, v1y);
		let deltaTheta = this.vectorAngle(v1x, v1y, v2x, v2y);

		if (fs === 0 && deltaTheta > 0) {
			deltaTheta -= Path.TAU;
		}
		if (fs === 1 && deltaTheta < 0) {
			deltaTheta += Path.TAU;
		}

		return [cx, cy, theta1, deltaTheta];
	}

	//
	// Approximate one unit arc segment with bézier curves,
	// see http://math.stackexchange.com/questions/873224
	//
	public static approximateUnitArc(theta1, deltaTheta) {
		const alpha = 4 / 3 * Math.tan(deltaTheta / 4);

		const x1 = Math.cos(theta1);
		const y1 = Math.sin(theta1);
		const x2 = Math.cos(theta1 + deltaTheta);
		const y2 = Math.sin(theta1 + deltaTheta);

		return [x1, y1, x1 - y1 * alpha, y1 + x1 * alpha, x2 + y2 * alpha, y2 - x2 * alpha, x2, y2];
	}

	public static a2c(x1, y1, x2, y2, fa, fs, rx, ry, phi) {
		const sinPhi = Math.sin(phi * Path.TAU / 360);
		const cosPhi = Math.cos(phi * Path.TAU / 360);

		// Make sure radii are valid
		//
		const x1p = cosPhi * (x1 - x2) / 2 + sinPhi * (y1 - y2) / 2;
		const y1p = -sinPhi * (x1 - x2) / 2 + cosPhi * (y1 - y2) / 2;

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
		const cc = Path.getArcCenter(x1, y1, x2, y2, fa, fs, rx, ry, sinPhi, cosPhi);

		const result = [];
		let theta1 = cc[2];
		let deltaTheta = cc[3];

		// Split an arc to multiple segments, so each segment
		// will be less than τ/4 (= 90°)
		//
		const segments = Math.max(Math.ceil(Math.abs(deltaTheta) / (Path.TAU / 4)), 1);
		deltaTheta /= segments;

		for (let i = 0; i < segments; i++) {
			result.push(Path.approximateUnitArc(theta1, deltaTheta));
			theta1 += deltaTheta;
		}

		// We have a bezier approximation of a unit circle,
		// now need to transform back to the original ellipse
		//
		return result.map((curve) => {
			for (let i = 0; i < curve.length; i += 2) {
				let x = curve[i + 0];
				let y = curve[i + 1];

				// scale
				x *= rx;
				y *= ry;

				// rotate
				const xp = cosPhi * x - sinPhi * y;
				const yp = sinPhi * x + cosPhi * y;

				// translate
				curve[i + 0] = xp + cc[0];
				curve[i + 1] = yp + cc[1];
			}

			return curve;
		});
	}
}


/*

call Path.getBoundingBox(path) from the outer
with d path
*/