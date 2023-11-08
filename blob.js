class Blob {
	constructor(centerRest) {
		this.radius = BLOB_RADIUS;
		this.centerRest = centerRest; // original location

		// CREATE PARTICLES:
		this.BP = []; //blob particles
		this.n = BLOB_PARTICLES;
        this.restArea = 0;
		let v0 = vec2(random(-100, 100), random(200, 220));
		for (let i = 0; i < this.n; i++) {
			let xi = this.radius * cos(i / this.n * TWO_PI) + centerRest.x;
			let yi = this.radius * sin(i / this.n * TWO_PI) + centerRest.y;
			let particle = createParticle(xi, yi);
			particle.v.set(v0);
			this.BP.push(particle);
		}
        for (let i = 0; i < this.n; i++) { 
			let p1 = this.BP[i];
			let p2 = this.BP[(i + 1) % this.n];

			let num = p1.p.x*p2.p.y - p1.p.y*p2.p.x;
            this.restArea += num/2;
		} 

		// CREATE EDGES FOR STRETCH SPRINGS + COLLISIONS:
		this.BE = []; // blob edges
		for (let i = 0; i < this.n; i++) {
			let p0 = this.BP[i];
			let p1 = this.BP[(i + 1) % this.n];
			this.BE.push(createEdge(p0, p1));
		}

		// SETUP APPEARANCE/FACIAL ELEMENTS: 
		// TODO
		let dc = 26;
		this.fillColor = color([221 + random(-dc, dc), 160 + random(-dc, dc), 221 + random(-dc, dc), 255]); // ("Plum"); // 221, 160, 221
	}

	blobParticles() {
		return this.BP;
	}

	// Loops over blob edges and accumulates stretch forces (Particle.f += ...)
	gatherForces_Stretch() {
		let k = STIFFNESS_STRETCH;
        let i = 0;
		for (let edge of this.BE) {
			let x0 = edge.lengthRest();
            let x = edge.length();
            let F = -k*(x-x0);
            
            let p1 = this.BP[i];
			let p2 = this.BP[(i + 1) % this.n];

            let forceX = F* (p1.p.x - p2.p.x)/ x; 
            let forceY = F* (p1.p.y - p2.p.y)/ x;
            
            p1.f.x += forceX;
            p1.f.y += forceY;

            p2.f.x -= forceX;
            p2.f.y -= forceY;
            i += 1;
		}
	}
	// Loops over blob particles and accumulates bending forces (Particle.f += ...)
	gatherForces_Bend() {
		let k = STIFFNESS_BEND;
		for (let i = 0; i < this.n; i++) { 
            
			let p0 = this.BP[(i-1 + this.n) % this.n];
			let p1 = this.BP[i];
			let p2 = this.BP[(i + 1) % this.n];
            
            let angle = TWO_PI/this.n;

			let a = createVector(p1.p.x - p0.p.x, p1.p.y - p0.p.y);
            let b = createVector(p2.p.x - p1.p.x, p2.p.y - p1.p.y);
            let a_hat = a.copy().normalize();
            let b_hat = b.copy().normalize();
            let F0Direction = b_hat.copy().sub(a_hat.copy().mult(a_hat.dot(b_hat)));
            let F0Magnitude = -k / (2 * a.mag());
            let F0 = F0Direction.mult(F0Magnitude);
            let F2Direction = a_hat.copy().sub(b_hat.copy().mult(a_hat.dot(b_hat)));
            let F2Magnitude = k / (2 * b.mag());
            let F2 = F2Direction.mult(F2Magnitude);
            let F1 = p5.Vector.add(F0, F2).mult(-1);
            
            p0.f.x += F0.x;
            p0.f.y += F0.y;
            p1.f.x += F1.x;
            p1.f.y += F1.y;
            p2.f.x += F2.x;
            p2.f.y += F2.y;
		}
	}
	// Loops over blob particles and gathers area compression forces (Particle.f += ...)
	gatherForces_Area() {
		let k = STIFFNESS_AREA;
        let A = 0;
		for (let i = 0; i < this.n; i++) { 
            
			let p_1 = this.BP[i];
			let p_2 = this.BP[(i + 1) % this.n];

			let num = p_1.p.x*p_2.p.y - p_1.p.y*p_2.p.x;
            A += num/2;
		} 
        for (let i = 0; i < this.n; i++) { 
            
            let p0 = this.BP[(i - 1 + this.n) % this.n];
			let p1 = this.BP[i];
			let p2 = this.BP[(i + 1) % this.n];
            let b = createVector(p2.p.x - p0.p.x, p2.p.y - p0.p.y);
            let grad = createVector(-b.y/2, b.x/2);
            
            let F = createVector(k*(A-this.restArea)*grad.x,k*(A-this.restArea)*grad.y);
            p1.f.x += F.x;
            p1.f.y += F.y;
		}
        
            
	}

	// Center of mass of all blob particles
	centerOfMass() {
		let com = vec2(0, 0);
		for (let particle of this.BP)
			acc(com, 1 / this.BP.length, particle.p); // assumes equal mass
		return com;
	}

	// Center of velocity of all blob particles
	centerOfVelocity() {
		let cov = vec2(0, 0);
		for (let particle of this.BP)
			acc(cov, 1 / this.BP.length, particle.v); // assumes equal mass
		return cov;
	}

	// Something simple to keep rigid blobs inside the box:
	rigidBounceOnWalls() {
		let pos = this.centerOfMass();
		let vel = this.centerOfVelocity();

		let R = BLOB_RADIUS + PARTICLE_RADIUS;

		// Boundary reflection (only if outside domain AND still moving outward):
		if ((pos.x < R && vel.x < 0) ||
			(pos.x > width - R && vel.x > 0)) {
			for (let particle of this.BP)
				particle.v.x *= -0.4;
		}
		if ((pos.y < R && vel.y < 0) ||
			(pos.y > height - R && vel.y > 0)) {
			for (let particle of this.BP)
				particle.v.y *= -0.4;
		}
	}

	// Something simple to keep nonrigid blob particles inside the box:
	nonrigidBounceOnWalls() {
		let R = PARTICLE_RADIUS;
		for (let particle of this.BP) {
			let pos = particle.p;
			let vel = particle.v;

			// Boundary reflection (only if outside domain AND still moving outward):
			if ((pos.x < R && vel.x < 0) ||
				(pos.x > width - R && vel.x > 0)) {
				vel.x *= -0.4;
			}
			if ((pos.y < R && vel.y < 0) ||
				(pos.y > height - R && vel.y > 0)) {
				vel.y *= -0.4;
			}
		}
	}

	draw() {
		push();
		strokeWeight(PARTICLE_RADIUS);
		stroke("DarkOrchid"); //BlueViolet");
		fill(this.fillColor); { // draw blob
			beginShape(TESS);
			for (let particle of this.BP)
				vertex(particle.p.x, particle.p.y);
			endShape(CLOSE);
		}

		if (DRAW_BLOB_PARTICLES) {
			fill("DarkOrchid");
			for (let particle of this.BP)
				circle(particle.p.x, particle.p.y, PARTICLE_RADIUS);
		}

		this.drawBlobFace();
		pop();
	}

	drawBlobFace() {
		push();
		// TODO: Draw your face here! :D 

		// CENTER OF MASS eyeball for now :/
		let com = this.centerOfMass();
		let cov = this.centerOfVelocity();
		stroke(0);
		fill(255);
		circle(com.x, com.y, 5 * PARTICLE_RADIUS);
		fill(0);
		circle(com.x + 0.01 * cov.x + 3 * sin(nTimesteps / 3), com.y + 0.01 * cov.y + random(-1, 1), PARTICLE_RADIUS);
		pop();
	}

}

// Creates edge and adds to edge list
function createEdge(particle0, particle1) {
	let edge = new Edge(particle0, particle1);
	edges.push(edge);
	return edge;
}

// Edge spring
class Edge {
	// Creates edge spring of default stiffness, STIFFNESS_STRETCH
	constructor(particle0, particle1) {
		this.q = particle0;
		this.r = particle1;
		this.restLength = this.q.pRest.dist(this.r.pRest);
		this.stiffness = STIFFNESS_STRETCH;
	}
	// True if both particles are pinned
	isRigid() {
		return (this.q.pin && this.r.pin);
	}
	// Current length of edge spring
	length() {
		return this.q.p.dist(this.r.p);
	}
	// Rest length of edge spring
	lengthRest() {
		return this.restLength;
	}
	// Draws the unstylized line 
	draw() {
		let a = this.q.p;
		let b = this.r.p;
		line(a.x, a.y, b.x, b.y);
	}
}
