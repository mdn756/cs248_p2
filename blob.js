class Blob {
	constructor(centerRest, ID) {
		this.radius = BLOB_RADIUS;
		this.centerRest = centerRest; // original location
		this.ID = ID;
		// CREATE PARTICLES:
		this.BP = []; //blob particles
		this.n = BLOB_PARTICLES;
        this.restArea = 0;
		let v0 = vec2(random(-100, 100), random(50, 70));
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
			let dampingForceX = -0.9 * (p1.v.x - p2.v.x);
			let dampingForceY = -0.9 * (p1.v.y - p2.v.y);
			forceX += dampingForceX;
			forceY += dampingForceY;	
            
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
            let grad = createVector(-b.y/2.0, b.x/2.0);
            
            let F = createVector(k*(A-this.restArea)*grad.x,k*(A-this.restArea)*grad.y);
            p1.f.x += F.x;
            p1.f.y += F.y;
		}
		//subtract off const times velocity
		// damping stretch spring
        
            
	}
	
	getAABB() {
		let minX = 999999;
		let minY = 999999;
		let maxX = 0;
		let maxY = 0;

		//with next point added
		for (let i = 0; i < this.n; i++) {  
			if (this.BP[i].p.x + deltaTime/1000*this.BP[i].v.x < minX) {
				minX = this.BP[i].p.x + deltaTime/1000*this.BP[i].v.x; 
			}
			if (this.BP[i].p.x < minX) {
				minX = this.BP[i].p.x ; 
			}
			if (this.BP[i].p.x + deltaTime/1000*this.BP[i].v.x > maxX) {
				maxX = this.BP[i].p.x + deltaTime/1000*this.BP[i].v.x; 
			}
			if (this.BP[i].p.x > maxX) {
				maxX = this.BP[i].p.x; 
			}
			if (this.BP[i].p.y + deltaTime/1000*this.BP[i].v.y < minY) {
				minY = this.BP[i].p.y + deltaTime/1000*this.BP[i].v.y; 
			}
			if (this.BP[i].p.y < minY) {
				minY = this.BP[i].p.y; 
			}
			if (this.BP[i].p.y + deltaTime/1000*this.BP[i].v.y > maxY) {
				maxY = this.BP[i].p.y + deltaTime/1000*this.BP[i].v.y; 
			}
			if (this.BP[i].p.y > maxY) {
				maxY = this.BP[i].p.y; 
			}
		} 
		
		let AABB = {
			minX: minX,
			minY: minY,
			maxX: maxX,
			maxY: maxY
		} 
		return AABB;
	}

	drawAABB() {
		/////////////////////TODO: COMMENT THIS IF YOU DONT WANT THE BV DRAWN 
		let AABB = this.getAABB();
		push();
		noFill();
		rect(AABB.minX, AABB.minY, AABB.maxX - AABB.minX, AABB.maxY - AABB.minY);
		pop();
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
		//this.drawAABB();
		pop();
	}

	drawBlobFace() {
		push();

		// CENTER OF MASS eyeball for now :/
		let com = this.centerOfMass();
		let cov = this.centerOfVelocity();
		if (this.ID == 0) {
			let x = com.x - seriousImage.width / 2;
    		let y = com.y - seriousImage.height / 2;
			let collisionDetected = false;
    		 
			let BV1 = blobs[0].getAABB();
			for (let i = 0; i < blobs.length; i++) {
				if (i !== this.ID) {
					let BV2 = blobs[i].getAABB();
					if (BV1.maxX >= BV2.minX && BV1.minX <= BV2.maxX) { //check x
						if (BV1.maxY >= BV2.minY && BV1.minY <= BV2.maxY) { //check y
							collisionDetected = true;
							
							break;
						}
					}
				}
			}
			if (collisionDetected) {
				image(yaranaika, x+WIDTH/75, y+HEIGHT/75);
			}
			else {
				image(seriousImage, x, y);
			}

			//TODO: on collision switch to yaranaika
		}
		else if (this.ID !== 0) {
			stroke(0);
			fill(255);
			circle(com.x-(BLOB_RADIUS*0.5), com.y, 5 * PARTICLE_RADIUS);
			fill(0);
			circle(com.x + 0.01 * cov.x + 3 * sin(nTimesteps / 3)-(BLOB_RADIUS*0.5), com.y + 0.01 * cov.y + random(-1, 1), PARTICLE_RADIUS);

			stroke(0);
			fill(255);
			circle(com.x+(BLOB_RADIUS*0.5), com.y, 5 * PARTICLE_RADIUS);
			fill(0);
			circle(com.x + 0.01 * cov.x + 3 * sin(nTimesteps / 3)+(BLOB_RADIUS*0.5), com.y + 0.01 * cov.y + random(-1, 1), PARTICLE_RADIUS);
		
			noFill();
			arc(com.x, com.y+(BLOB_RADIUS*0.25), 5 * PARTICLE_RADIUS, 5 * PARTICLE_RADIUS, 0, PI);
			pop();
		}
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
		
		//let AABB = this.getBV()
		//noFill();
		//stroke(0, 0, 255);
		//rect(AABB.minX, AABB.minY, AABB.maxX - AABB.minX, AABB.maxY - AABB.minY);

	}
	getBV() {
		let a = this.q.p;
		let b = this.r.p;
		let d = 10;
		let minX = min(a.x, b.x) - d;
		let minY = min(a.y, b.y) - d;
		let maxX = max(a.x, b.x) + d;
		let maxY = max(a.y, b.y) + d;
		let AABB = {
			minX: minX,
			minY: minY,
			maxX: maxX,
			maxY: maxY
		} 
		return AABB;
	}
}
