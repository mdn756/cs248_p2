class SDF {

	/**
	 * Computes SDF union distance value, min(d1,d2), and its corresponding object reference.
	 * @param distO1 {[number,object]} Obstacle #1's SDF value, d, and object ref, o, as an array [d1,o1].
	 * @param distO2 {[number,object]} Obstacle #2's SDF value, d, and object ref, o, as an array [d2,o2].
	 * @return {[number,object]} Returns [d1,o1] if d1<d2, and [d2,o2] otherwise.
	 */
	static opUnion(distO1, distO2) {
		return (distO1[0] < distO2[0]) ? distO1 : distO2;
	}
	/**
	 * Computes SDF intersection distance value, max(d1,d2), and its corresponding object reference.
	 * @param distO1 {[number,object]} Obstacle #1's SDF value, d, and object ref, o, as an array [d1,o1].
	 * @param distO2 {[number,object]} Obstacle #2's SDF value, d, and object ref, o, as an array [d2,o2].
	 * @return {[number,object]} Returns [d1,o1] if d1>d2, and [d2,o2] otherwise.
	 */
	static opIntersection(distO1, distO2) {
		return (distO1[0] > distO2[0]) ? distO1 : distO2;
	}

	// ANOTHER BOOLEAN YOU MAY WANT TO USE:
	//float opSubtraction ( float d1, float d2 ) { return max(-d1,d2);}

	/**
	 * SDF of a circle at origin.
	 * @param {p5.Vector} p Evaluation point
	 * @param {number}    r Radius
	 * @return {number}   SDF value
	 */
	static sdCircle(p, r) {
		return length(p) - r;
	}

	/**
	 * SDF of a box at origin
	 * @param {p5.Vector} p Evaluation point
	 * @param {p5.Vector} b Half widths in X & Y
	 * @return {number}   SDF value
	 */
	static sdBox(p, b) {
		let d = sub(absv2(p), b);
		return length(maxv2(d, 0.0)) + min(max(d.x, d.y), 0.0);
	} 

	static sdRBox(p, a, b, r) {
		let ba = createVector(b.x - a.x, b.y - a.y);
		let pa = createVector(p.x - a.x, p.y - a.y);
		let h = constrain(pa.dot(ba) / ba.magSq(), 0, 1);
		let q = p5.Vector.sub(pa, ba.mult(h));
		let d = q.mag() - r;  // Subtract the radius from the distance
		return d;
	}

	
	static sdArc(p) {
        let s = 0.9;
        let c = 0.44;
        let sc = createVector(s,c);
        p.x = abs(p.x);
        if (c*p.x>s*p.y){
            return length(sub(p,mult(sc, 1)))
        } 
        else {
            return abs(length(p)-1)-0.1;
        }
	}
	
	
	
	static sdTriangle(p, p0, p1, p2) {
		let e0 = sub(p1, p0);
		let e1 = sub(p2, p1);
		let e2 = sub(p0, p2);

		let v0 = sub(p, p0);
		let v1 = sub(p, p1);
		let v2 = sub(p, p2);

		let pq0 = sub(v0, mult(e0,clamp( dot(v0,e0)/dot(e0,e0), 0.0, 1.0 )));
		let pq1 = sub(v1, mult(e1,clamp( dot(v1,e1)/dot(e1,e1), 0.0, 1.0 )));
		let pq2 = sub(v2, mult(e2,clamp( dot(v2,e2)/dot(e2,e2), 0.0, 1.0 )));

		let s = sign(e0.x*e2.y - e0.y*e2.x);
		let s0 = createVector( dot( pq0, pq0 ), s*(v0.x*e0.y-v0.y*e0.x) );
		let s1 = createVector( dot( pq1, pq1 ), s*(v1.x*e1.y-v1.y*e1.x) );
		let s2 =  createVector( dot( pq2, pq2 ), s*(v2.x*e2.y-v2.y*e2.x) );
		let xarr = [s0.x, s1.x, s2.x];
		let yarr = [s0.y, s1.y, s2.y];
		let dx = min(xarr);
		let dy = min(yarr);
		return -1 * Math.sqrt(dx)*Math.sign(dy);
	}
	
	// static sdRoundBox(p, b) {
	// 	let d = sub(absv2(p), b);
	// 	return length(maxv2(d, 0.0)) + min(max(d.x, d.y), 0.0);
	// }

	/** 
	 * Uneven Capsule - exact   (https://www.shadertoy.com/view/4lcBWn)
	 * @param {p5.Vector} q  Evaluation point
	 * @param {number}    r1 
	 * @param {number}    r2
	 * @param {number}    h  
	 */
	static sdUnevenCapsule(q, r1, r2, h) {
		// p.x = abs(p.x);
		let p = vec2(abs(q.x), q.y); // local copy to avoid changing pass-by-ref input
		const b = (r1 - r2) / h;
		const a = Math.sqrt(1.0 - b * b);
		const k = dot(p, vec2(-b, a));
		if (k < 0.0) return length(p) - r1;
		if (k > a * h) return length(sub(p, vec2(0.0, h))) - r2;
		return dot(p, vec2(a, b)) - r1;
	}
}