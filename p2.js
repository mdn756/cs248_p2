/** 
 * StarterCode for "Attack of the Blobs!" 
 * CS248B Fundamentals of Computer Graphics: Animation & Simulation
 * 
 * Fill in the the missing code (see TODO items).
 * Try reducing MAX_BLOBS to 1 to get started. 
 * Good luck!!
 * 
 * @author Doug L. James <djames@cs.stanford.edu> 
 * Sarah Teaw sarahlst@stanford.edu
 * Maxime Nee mdn756@stanford.edu
 * @date 10/28/2022
 */

const MAX_BLOBS = 3; /// TODO: 100 or more to complete "Attack of the Blobs!" challenge. Use just a few for testing. 
const DRAW_BLOB_PARTICLES = true;

const WIDTH = 1024;
const HEIGHT = 1024;
const PARTICLE_RADIUS = WIDTH / 400.0; // for rendering
const PARTICLE_MASS = 1.0;
const BLOB_PARTICLES = 12; // 12 (F22)
const BLOB_RADIUS = WIDTH / 25 ;

const STIFFNESS_STRETCH = 10000.0; // TODO: Set as you wish
const STIFFNESS_BEND = 100000.0; //    TODO: Set as you wish
const STIFFNESS_AREA = 1 * 3.14*BLOB_RADIUS*BLOB_RADIUS/5000; //    TODO: Set as you wish
//////// IMPORTANT ARRAYS OF THINGS /////////
let particles = []; // All particles in the scene (rigid + blobs) 
let edges = []; //     All edges in the scene (rigid + blobs)
let blobs = []; //     All blobs in the scene (increases over time)
let environment; //    Environment with all rigid edges available as getEdges()

let isPaused = true;
let nTimesteps = 0; // #frame-length timesteps taken
let detectedEdgeEdgeFailure = false; // Halts simulation and turns purple if true -- blobs win!
let ID = 0;
// Graph paper texture map:
let bgImage;

function preload() {
	bgImage = loadImage('graphpaper.jpg');
	seriousImage = loadImage('serious.png');
	yaranaika = loadImage('yaranaika.png');
}

function setup() {
	createCanvas(WIDTH, HEIGHT);
	seriousImage.resize(int(seriousImage.width*8 / BLOB_RADIUS), int(seriousImage.height*8 / BLOB_RADIUS));
	yaranaika.resize(int(yaranaika.width*6 / BLOB_RADIUS), int(yaranaika.height*6 / BLOB_RADIUS));
	background(100);
	ellipseMode(RADIUS);
	environment = new Environment();
	//print("|particles|=" + particles.length + ",  |edge|=" + edges.length + ",  |blobs|=" + blobs.length);
}

/// Timesteps (w/ substeps) and draws everything.
function draw() {

	///// SIMULATE /////
	if (!isPaused) {
		// CREATE BLOBS 
		if (nTimesteps % 10 == 0) { 
			if (blobs.length < MAX_BLOBS) 
				createRandomBlob(ID); // tries to create one if free space available
				ID += 1;
		}

		// TIMESTEP!
		let dtFrame = 0.01;
		let nSubsteps = 1; // #times to split dtFrame
		for (let step = 0; step < nSubsteps; step++)
			advanceTime(dtFrame / nSubsteps);
		nTimesteps++;
	}

	///// RENDER /////
	push();
	background(0);
	environment.draw();
	for (let blob of blobs)
		blob.draw();
	pop();
	drawMouseForce();

	/// TEXT OUTPUT:
	push();
	textSize(18);
	noStroke();
	fill(0);
	text("#BLOBS: " + blobs.length, 10, 20);
	text("#EDGES: " + edges.length, 10, 40);
	text("#PARTICLES: " + particles.length, 10, 60);
	pop();
}

function keyPressed() {
	if (keyCode == 32) // spacebar
		isPaused = !isPaused;
}

function advanceTime(dt) {
	environment.advanceTime(dt);

	//////////////////////////////////////
	////// GATHER PARTICLE FORCES ////////
	{
		// Clear forces:
		for (let particle of particles)
			particle.f.set(0, 0);

		gatherParticleForces_Gravity();

		// Damping (springs or otherwise -- you can add some if you want): 
                                                         
		// Blob springs: 
		for (let blob of blobs) {
			blob.gatherForces_Stretch();
			blob.gatherForces_Bend();
			blob.gatherForces_Area();
		}


		gatherParticleForces_Penalty();

		// Mouse force (modify if you want):
		applyMouseForce();
	}

	//////////////////////////////////////////
	// Update velocity (using mass filtering):
	for (let particle of particles)
		acc(particle.v, dt * particle.invMass(), particle.f)

	//////////////////////////////////////////
	// Collision filter: Correct velocities //
	applyPointEdgeCollisionFilter();

	//////////////////////////////////////////
	// Update positions:
	for (let particle of particles)
		acc(particle.p, dt, particle.v)

	verifyNoEdgeEdgeOverlap();
	isBVOverlap();
}

function isBVOverlap() {
	//blob against blob
	for (let i = 0; i < blobs.length; i++) {
		for (let j = i+1; j < blobs.length; j++) {
			BV1 = blobs[i].getAABB();
			BV2 = blobs[j].getAABB();
			if (BV1.maxX >= BV2.minX && BV1.minX <= BV2.maxX) { //check x
				if (BV1.maxY >= BV2.minY && BV1.minY <= BV2.maxY) { //check y
					//isPaused = true;
					return true;
				}
			}
		}
	}
	// blob against env edges
	for (let i = 0; i < blobs.length; i++) {
		for (let j = 1; j < environment.envEdges.length; j++) {
			BV1 = blobs[i].getAABB();
			BV2 = environment.envEdges[j].getBV();
			if (BV1.maxX >= BV2.minX && BV1.minX <= BV2.maxX) { //check x
				if (BV1.maxY >= BV2.minY && BV1.minY <= BV2.maxY) { //check y
					//isPaused = true;
					return true;
				}
			}
		}
	}

	return false;
}

function isCollision(t, particle, r, q) {
	let delta_t = 0.01;
	if (t < delta_t) {
		let pt = particle.p.copy();
		let rt = r.p.copy();
		let qt = q.p.copy();
		acc(pt, t, particle.v);
		acc(rt, t, r.v);
		acc(qt, t, q.v);
		let onX = (pt.x >= (rt.x - 2.0) && pt.x <= (qt.x+2.0)) || (pt.x >= (qt.x-2.0) && pt.x <= (rt.x+2.0));
		let onY = (pt.y >= (rt.y - 2.0) && pt.y <= (qt.y+2.0)) || (pt.y >= (qt.y-2.0) && pt.y <= (rt.y+2.0));
		if (onX && onY) {
			return true;
		}
	}
	return false;
}

function applyPointEdgeCollisionFilter() {
	let epsilon = 0.5;
	// TEMP HACK (remove!): rigid bounce off walls so they don't fly away
	for (let blob of blobs) blob.nonrigidBounceOnWalls();

	// TODO: Process all point-edge CCD impulses 
	// FIRST: Just rigid edges.
	for (let edge of environment.getEdges()) {
		// colinearity check to find a time and position
		for (let blob of blobs) {
			for (let particle of blob.BP){
				let q = edge.q;
				let r = edge.r;
				let A = sub(r.p, q.p);
				let B = sub(particle.p, q.p);
				let a_dot = sub(r.v, q.v);
				let b_dot = sub(particle.v, q.v);
				let a_dot_w_b_dot = (a_dot.x * b_dot.y) - (a_dot.y * b_dot.x);
				let a_dot_w_b = (a_dot.x * B.y) - (a_dot.y * B.x);
				let a_w_b_dot = (A.x * b_dot.y) - (A.y * b_dot.x);
				let a_w_b = (A.x * B.y) - (A.y * B.x);
				let a = a_dot_w_b_dot;
				let b = a_dot_w_b + a_w_b_dot;
				let c = a_w_b;
				let discriminant = sq(b) - (4.0 * a * c);
				let t = 0;
				if (a == 0) {
					t = -1.0*c/b;
					let delta_t = 0.01;
					if (t < delta_t && t > 0) {
						let pt = particle.p.copy();
						let rt = r.p.copy();
						let qt = q.p.copy();
						acc(pt, t, particle.v);
						acc(rt, t, r.v);
						acc(qt, t, q.v);
						let onX = (pt.x >= (rt.x - 2.0) && pt.x <= (qt.x+2.0)) || (pt.x >= (qt.x-2.0) && pt.x <= (rt.x+2.0));
						let onY = (pt.y >= (rt.y - 2.0) && pt.y <= (qt.y+2.0)) || (pt.y >= (qt.y-2.0) && pt.y <= (rt.y+2.0));
						if (onX && onY) {
							let pr = sub(pt, rt).mag();
							let pq = sub(qt, pt).mag();
							let rq = sub(rt, qt).mag();
							let alpha = pr/rq;
							let beta = pq/rq;
							let c_dot = q.p.copy();
							acc(c_dot, sub(r.v, q.v), alpha);
							let contact_v = sub(pt, c_dot);
							let n0 = createVector(q.p.y-r.p.y, r.p.x-q.p.x);
							n0.normalize();
							let n_hat = mult(n0, sign(-1*dot(contact_v, n0)));
							let v_n = dot(sub(particle.v, c_dot), n_hat);
							let m_eff = 1.0/particle.mass;
							let gamma = v_n*-1.0*(1+epsilon)*m_eff;
							acc(particle.v, gamma/particle.mass, n_hat);
						}
					}
				}
			}
		}
	}
	// SECOND: All rigid + blob edges (once you get this ^^ working)
	// edgesToCheck = edges;
	for (let blob of blobs) {
		// colinearity check to find a time and position
		for (let edge of blob.BE) {
			for (let particle of environment.getParticles()){
				let q = edge.q;
				let r = edge.r;
				let A = sub(r.p, q.p);
				let B = sub(particle.p, q.p);
				let a_dot = sub(r.v, q.v);
				let b_dot = sub(particle.v, q.v);
				let a_dot_w_b_dot = (a_dot.x * b_dot.y) - (a_dot.y * b_dot.x);
				let a_dot_w_b = (a_dot.x * B.y) - (a_dot.y * B.x);
				let a_w_b_dot = (A.x * b_dot.y) - (A.y * b_dot.x);
				let a_w_b = (A.x * B.y) - (A.y * B.x);
				let a = a_dot_w_b_dot;
				let b = a_dot_w_b + a_w_b_dot;
				let c = a_w_b;
				let discriminant = sq(b) - (4.0 * a * c);
				let t = 0;
				let delta_t = 0.01;
				if (a == 0) {
					t = -1.0*c/b;
				}
				if (b == 0) {
					if (c > 0) return;
					t = sqrt(-c/a);
				}
				if (discriminant > 0) {
					let factor = -0.5 * (b + sign(b)*sqrt(discriminant));
					let t1 = factor/a;
					let t2 = c/factor;
					if (isCollision(t1, particle, r, q)) t = t1;
					else if (isCollision(t2, particle, r, q)) t = t2;
					else return;
				}
				if (t < delta_t && t > 0) {
					let pt = particle.p.copy();
					let rt = r.p.copy();
					let qt = q.p.copy();
					acc(pt, t, particle.v);
					acc(rt, t, r.v);
					acc(qt, t, q.v);
					let onX = (pt.x >= (rt.x - 2.0) && pt.x <= (qt.x+2.0)) || (pt.x >= (qt.x-2.0) && pt.x <= (rt.x+2.0));
					let onY = (pt.y >= (rt.y - 2.0) && pt.y <= (qt.y+2.0)) || (pt.y >= (qt.y-2.0) && pt.y <= (rt.y+2.0));
					if (onX && onY) {
						let pr = sub(pt, rt).mag();
						let pq = sub(qt, pt).mag();
						let rq = sub(rt, qt).mag();
						let alpha = pr/rq;
						let beta = pq/rq;
						let c_dot = q.p.copy();
						acc(c_dot, sub(r.v, q.v), alpha);
						let contact_v = sub(pt, c_dot);
						let n0 = createVector(q.p.y-r.p.y, r.p.x-q.p.x);
						n0.normalize();
						let n_hat = mult(n0, sign(-1*dot(contact_v, n0)));
						let v_n = dot(sub(particle.v, c_dot), n_hat);
						let m_eff = 0;
						if (!particle.pin) m_eff += 1.0/particle.mass;
						if (!r.pin) m_eff += beta*beta/r.mass;
						if (!q.pin) m_eff += alpha*alpha/q.mass;
						let gamma = v_n*-1.0*(1+epsilon)*(m_eff);
						acc(q.v, -beta*gamma/q.mass, n_hat);
						acc(r.v, -alpha*gamma/r.mass, n_hat);
					}
				}	
			}

		}
	}

	// Initially just brute force all-pairs checks, later use bounding volumes or better broad phase.

}

// Efficiently checks that no pair of edges overlap, where the pairs do not share a particle in common.
function verifyNoEdgeEdgeOverlap() {
	if (detectedEdgeEdgeFailure) return; // already done

	// TODO: Optional: Make faster with broad phase
	// SIMPLE: Brute force check on edges i<j:
	for (let i = 0; i < edges.length - 1; i++) {
		let ei = edges[i];
		for (let j = i + 1; j < edges.length; j++) {
			let ej = edges[j];
			if (checkEdgeEdgeOverlap(ei, ej)) {
				// HALT!
				detectedEdgeEdgeFailure = true;
				isPaused = true;
				return;
			}
		}
	}
}

function checkEdgeEdgeOverlap(ei, ej) {
	// TODO: Implement robust-enough test. 
	let a = ei.q.p;
	let b = ei.r.p;
	let c = ej.q.p;
	let d = ej.r.p;

	if (triSign(c,d,a)*triSign(c,d,b) < 0) {
		if (triSign(a,b,c)*triSign(a,b,d) < 0) {
			return true;
		}
	}
	return false;
}

function triSign(a, b, c) {
	let dx1 = b.x - a.x;
	let dx2 = c.x - b.x;
	let dy1 = b.y - a.y;
	let dy2 = c.y - b.y;

	let sign = -dy1*dx2 + dy2*dx1;
	return sign
}

// Computes penalty forces between all point-edge pairs
function gatherParticleForces_Penalty() {
	let k =1000.0;
	// let warmup = true;
	// if (warmup) { // First just consider rigid environment edges:
	for (let blob of blobs) {
		for (let edge of environment.getEdges()) {
			// TODO (part1): Apply point-edge force (if pt not on edge!)
			for (let particle of blob.BP) {
				// Calculate the vector from the particle to the edge
				let q = edge.q;
				let r = edge.r;
				let edgeVector = sub(r.p, q.p);
				let particleToEdge = sub(particle.p, q.p);
				
				// Calculate the projection of particleToEdge onto edgeVector
				let projection = edgeVector.dot(particleToEdge) / edgeVector.magSq();
				
				// Calculate the closest point on the edge to the particle
				let closestPoint;
				projection = constrain(projection, 0, 1);
				closestPoint = add(q.p, p5.Vector.mult(edgeVector, projection));
				
				// Calculate the distance between the particle and the closest point on the edge
				let n = sub(particle.p, closestPoint);
				let distance = n.mag();

				// Specify a threshold distance, below which penalty forces are applied
				let d0 = 15.0; // Adjust this threshold as needed
				
				// Apply a penalty force if the particle is too close to the edge
				if (distance < d0) {
					let penaltyForceMagnitude = (d0 - distance) * (k); // Adjust 'k' as needed
					let penaltyForce = createVector(particle.p.x-closestPoint.x, particle.p.y-closestPoint.y);
					penaltyForce.normalize();
					penaltyForce.mult(penaltyForceMagnitude);
					particle.f.add(penaltyForce);
				}
            }
		}
		for (let edge of blob.BE) {
			for (let particle of environment.getParticles()) {
				// Calculate the vector from the particle to the edge
				let q = edge.q;
				let r = edge.r;
				let edgeVector = sub(r.p, q.p);
				let particleToEdge = sub(particle.p, q.p);
				
				// Calculate the projection of particleToEdge onto edgeVector
				let projection = edgeVector.dot(particleToEdge) / edgeVector.magSq();
				
				// Calculate the closest point on the edge to the particle
				let closestPoint;
				projection = constrain(projection, 0, 1);
				closestPoint = add(q.p, p5.Vector.mult(edgeVector, projection));
				
				// Calculate the distance between the particle and the closest point on the edge
				let n = sub(particle.p, closestPoint);
				let distance = n.mag();

				// Specify a threshold distance, below which penalty forces are applied
				let d0 = 15.0; // Adjust this threshold as needed
				
				// Apply a penalty force if the particle is too close to the edge
				if (distance < d0) {
					let penaltyForceMagnitude = (d0 - distance) * (-k); // Adjust 'k' as needed
					let penaltyForce = createVector(particle.p.x-closestPoint.x, particle.p.y-closestPoint.y);
					penaltyForce.normalize();
					penaltyForce.mult(penaltyForceMagnitude);
					edge.q.f.add(penaltyForce);
					edge.r.f.add(penaltyForce);
				}
			}
		}
		for (let blob2 of blobs) {
			if (blob2 == blob) {
				continue;
			}
			for (let edge of blob2.BE) {
				// TODO (part1): Apply point-edge force (if pt not on edge!)
				for (let particle of blob.BP) {
					// Calculate the vector from the particle to the edge
					let q = edge.q;
					let r = edge.r;
					let edgeVector = sub(r.p, q.p);
					let particleToEdge = sub(particle.p, q.p);
					
					// Calculate the projection of particleToEdge onto edgeVector
					let projection = edgeVector.dot(particleToEdge) / edgeVector.magSq();
					
					// Calculate the closest point on the edge to the particle
					let closestPoint;
					projection = constrain(projection, 0, 1);
					closestPoint = add(q.p, p5.Vector.mult(edgeVector, projection));
					
					// Calculate the distance between the particle and the closest point on the edge
					let n = sub(particle.p, closestPoint);
					let distance = n.mag();
	
					// Specify a threshold distance, below which penalty forces are applied
					let d0 = 15.0; // Adjust this threshold as needed
					
					// Apply a penalty force if the particle is too close to the edge
					if (distance < d0) {
						let penaltyForceMagnitude = (d0 - distance) * (k); // Adjust 'k' as needed
						let penaltyForce = createVector(particle.p.x-closestPoint.x, particle.p.y-closestPoint.y);
						penaltyForce.normalize();
						penaltyForce.mult(penaltyForceMagnitude);
						particle.f.add(penaltyForce);
					}
				}
			}
	

		}
			// TODO (part2): Apply point-edge force (if pt not on edge!)
	}
	// }
}

function gatherParticleForces_Gravity() {
	let g = vec2(0, 100); //grav accel
	for (let particle of particles) {
		acc(particle.f, particle.mass, g); // f += m g
		// acc(particle.f, particle.v, -10.0);
	}
}



// Creates a default particle and adds it to particles list
function createParticle(x, y) {
	let p = new Particle(vec2(x, y), 1.0, PARTICLE_RADIUS);
	particles.push(p);
	return p;
}

class Particle {
	constructor(pRest, mass, radius) {
		this.pRest = vec2(pRest.x, pRest.y);
		this.p = vec2(pRest.x, pRest.y);
		this.v = vec2(0, 0);
		this.pin = false; // true if doesn't respond to forces
		this.mass = mass;
		this.radius = radius;
		this.f = vec2(0, 0);
	}
	invMass() {
		return (this.pin ? 0.0 : 1.0 / this.mass);
	}
	// Emits a circle
	draw() {
		// nobody = (this.pin ? fill("red") : fill(0)); // default colors (red if pinned)
		circle(this.p.x, this.p.y, this.radius); //ellipseMode(RADIUS);
	}
}

function drawAABB(ei) {
	// TODO: Implement robust-enough test. 
	let a = ei.q.p;
	let b = ei.r.p;

	return false;
}

// RIGID ENVIRONMENT COMPOSED OF LINE SEGMENTS (pinned Edges)
class Environment {

	constructor() {
		this.envParticles = [];
		this.envEdges = [];

		///// BOX /////
		let r = PARTICLE_RADIUS;
		this.p00 = createParticle(r, r);
		this.p01 = createParticle(r, HEIGHT - r);
		this.p11 = createParticle(WIDTH - r, HEIGHT - r);
		this.p10 = createParticle(WIDTH - r, r);
		this.p00.pin = this.p01.pin = this.p11.pin = this.p10.pin = true;
		this.envParticles.push(this.p00);
		this.envParticles.push(this.p01);
		this.envParticles.push(this.p11);
		this.envParticles.push(this.p10);
		this.envEdges.push(createEdge(this.p00, this.p01));
		this.envEdges.push(createEdge(this.p01, this.p11));
		this.envEdges.push(createEdge(this.p11, this.p10));
		this.envEdges.push(createEdge(this.p10, this.p00));

		///// OBSTACLES FOR FUN /////
		{
			// (F22) ANGLED LINES: 
			//for (let i = 0.5; i < 4; i++) this.createEnvEdge(i * width / 5, height / 2, (i + 1) * width / 5, height * 0.75);

			// (F23) PACHINKO PEGS:
			let n = 8;
			for (let i = 1; i < n; i++) {
				for (let j = 1; j < n; j++) {
					if ((i + j) % 2 == 0) { // alternating pegs
						//this.createPachinkoPeg(width * (i / n), height / 5 + height / 4 * 3 * (j / n), 11); // round 4-edge peg
						this.createPachinkoWedge(width * (i / n), height / 5 + height / 4 * 3 * (j / n), 11); // cheap 2-edge peg
					}
				}
			}
		}
	}

	// Returns all rigid-environment edges.
	getEdges() {
		return this.envEdges;
	}
	// Returns all rigid-environment edges.
	getParticles() {
		return this.envParticles;
	}
	// Creates a lone rigid edge.
	createEnvEdge(x0, y0, x1, y1) {
		let p0 = createParticle(x0, y0);
		let p1 = createParticle(x1, y1);
		p0.pin = true;
		p1.pin = true;
		let e = createEdge(p0, p1);
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envEdges.push(e);
	}
	// Create a lone roundish peg at (x0,y0) with radius r.
	createPachinkoPeg(x, y, r) {
		let p0 = createParticle(x + r, y);
		let p1 = createParticle(x, y + r);
		let p2 = createParticle(x - r, y);
		let p3 = createParticle(x, y - r);
		p0.pin = p1.pin = p2.pin = p3.pin = true;
		let e01 = createEdge(p0, p1);
		let e12 = createEdge(p1, p2);
		let e23 = createEdge(p2, p3);
		let e30 = createEdge(p3, p0);
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envParticles.push(p2);
		this.envParticles.push(p3);
		this.envEdges.push(e01);
		this.envEdges.push(e12);
		this.envEdges.push(e23);
		this.envEdges.push(e30);
	}
	// Create a lighter-weight wedge-shaped peg at (x0,y0) with radius r.
	createPachinkoWedge(x, y, r) {
		let p0 = createParticle(x + r, y);
		let p1 = createParticle(x, y - r);
		let p2 = createParticle(x - r, y);
		p0.pin = p1.pin = p2.pin = true;
		let e01 = createEdge(p0, p1);
		let e12 = createEdge(p1, p2);
		this.envParticles.push(p0);
		this.envParticles.push(p1);
		this.envParticles.push(p2);
		this.envEdges.push(e01);
		this.envEdges.push(e12);
			
	}

	// Updates any moveable rigid elements
	advanceTime(dt) {}

	// Makes popcorn <jk> no it doesn't... 
	draw() {
		push();
		image(bgImage, 0, 0, WIDTH, HEIGHT);

		if (detectedEdgeEdgeFailure) { // HALT ON OVERLAP + DRAW PURPLE SCREEN
			push();
			fill(191, 64, 191, 150);
			rect(0, 0, width, height);
			pop();
		}

		stroke("blue");
		strokeWeight(PARTICLE_RADIUS);
		for (let edge of this.envEdges) {
			edge.draw();
		}
		fill("blue");
		noStroke();
		for (let particle of this.envParticles) {
			particle.draw();
		}
		pop(); // wait, it does pop :/ 
	}
}

// Creates a blob centered at (x,y), and adds things to lists (blobs, edges, particles).
function createBlob(x, y, ID) {
	let b = new Blob(vec2(x, y), ID);
	blobs.push(b);
	return b;
}

// Tries to create a new blob at the top of the screen. 
function createRandomBlob(ID) {
	for (let attempt = 0; attempt < 5; attempt++) {
		let center = vec2(random(2 * BLOB_RADIUS, WIDTH - 2 * BLOB_RADIUS), BLOB_RADIUS * 2.3); //random horizontal spot
		// CHECK TO SEE IF NO BLOBS NEARBY:
		let tooClose = false;
		for (let blob of blobs) {
			let com = blob.centerOfMass();
			if (com.dist(center) < 3 * blob.radius) // too close
				tooClose = true;
		}
		// if we got here, then center is safe:
		if (!tooClose) {
			createBlob(center.x, center.y, ID);
			return;
		}
	}
}




/////////////////////////////////////////////////////////////////
// Some convenient GLSL-like macros for p5.Vector calculations //
/////////////////////////////////////////////////////////////////
function length(v) {
	return v.mag();
}

function dot(x, y) {
	return x.dot(y);
}

function dot2(x) {
	return x.dot(x);
}

function vec2(a, b) {
	return createVector(a, b);
}

function vec3(a, b, c) {
	return createVector(a, b, c);
}

function sign(n) {
	return Math.sign(n);
}

function clamp(n, low, high) {
	return constrain(n, low, high);
}

function add(v, w) {
	return p5.Vector.add(v, w);
}
function mult(v, w) {
	return p5.Vector.mult(v, w);
}
function sub(v, w) {
	return p5.Vector.sub(v, w);
}

function absv2(v) {
	return vec2(Math.abs(v.x), Math.abs(v.y));
}

function maxv2(v, n) {
	return vec2(Math.max(v.x, n), Math.max(v.y, n));
}

function minv2(v, n) {
	return vec2(Math.min(v.x, n), Math.min(v.y, n));
}

function vertexv2(p) {
	vertex(p.x, p.y);
}

// v += a*w
function acc(v, a, w) {
	v.x += a * w.x;
	v.y += a * w.y;
}

function rotateVec2(v, thetaRad) {
	const c = cos(thetaRad);
	const s = sin(thetaRad);
	return vec2(c * v.x - s * v.y, s * v.x + c * v.y);
}