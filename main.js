// Initialise the canvas
const canvas = document.querySelector('canvas');
const ctx = canvas.getContext('2d');


/** Helper function to draw a circle on the canvas */
function drawCircle(x, y, radius, color) {
    ctx.fillStyle = color;
    ctx.beginPath();
    ctx.arc(x, y, radius, 0, 2 * Math.PI);
    ctx.fill();
}

/** Helper function to draw a line between two points */
function drawLine(coords1, coords2) {
    ctx.strokeStyle = 'white';
    ctx.beginPath();
    ctx.moveTo(coords1.x, coords1.y);
    ctx.lineTo(coords2.x, coords2.y);
    ctx.stroke();
}

/** Calculates the equation of motion for the first pendulum */
function accel1() {
    return (
        (
            (-G * (2 * mass1.m + mass2.m) * Math.sin(mass1.the)) +
            (-mass2.m * G * Math.sin(mass1.the - 2 * mass2.the)) +
            (-2 * Math.sin(mass1.the - mass2.the) * mass2.m) *
            (mass2.v * mass2.v * length2 + mass1.v * mass1.v * length1 * Math.cos(mass1.the - mass2.the))
        ) / 
        (length1 * (2 * mass1.m + mass2.m - mass2.m * Math.cos(2 * mass1.the - 2 * mass2.the)))
    );
}

/** Calculates the equation of motion for the second pendulum */
function accel2() {
    return (
        (2 * Math.sin(mass1.the - mass2.the)) *
        (
            (mass1.v * mass1.v * length1 * (mass1.m + mass2.m)) +
            (G * (mass1.m + mass2.m) * Math.cos(mass1.the)) +
            (mass2.v * mass2.v * length2 * mass2.m * Math.cos(mass1.the - mass2.the))
        ) /
        (length2 * (2 * mass1.m + mass2.m - mass2.m * Math.cos(2 * mass1.the - 2 * mass2.the)))
    );
}


// Set the initial condition for the simulation
const MASS = 60;
const G = 9.81;
const length1 = 200;
const length2 = 200;

let pivot = {
    x: canvas.width / 2,
    y: canvas.height / 3
}

let mass1 = {
    m: MASS,
    the: Math.PI / 2,
    v: 0,
    a: 0
};

let coords1 = {
    x: pivot.x + length1 * Math.sin(mass1.the),
    y: pivot.y + length1 * Math.cos(mass1.the)
};

let mass2 = {
    m: MASS,
    the: Math.PI / 2,
    v: 0,
    a: 0
};

let coords2 = {
    x: coords1.x + length2 * Math.sin(mass2.the),
    y: coords1.y + length2 * Math.cos(mass2.the)
};

/** The main function which uses the Runge-Kutta method (described in the Python solution)
 * to update the position of the pendulums.
 * I won't go into the details of what each step does in this function, because it is the same as the Python solution.
 */
function update(dt) {

    let k1_v1 = dt * accel1(mass1.the, mass2.the, mass1.v, mass2.v);
    let k1_v2 = dt * accel2(mass1.the, mass2.the, mass1.v, mass2.v);
    let k1_the1 = dt * mass1.v;
    let k1_the2 = dt * mass2.v;

    let k2_v1 = dt * accel1(mass1.the + k1_the1 / 2, mass2.the + k1_the2 / 2, mass1.v + k1_v1 / 2, mass2.v + k1_v2 / 2);
    let k2_v2 = dt * accel2(mass1.the + k1_the1 / 2, mass2.the + k1_the2 / 2, mass1.v + k1_v1 / 2, mass2.v + k1_v2 / 2);
    let k2_the1 = dt * (mass1.v + k1_v1 / 2);
    let k2_the2 = dt * (mass2.v + k1_v2 / 2);

    let k3_v1 = dt * accel1(mass1.the + k2_the1 / 2, mass2.the + k2_the2 / 2, mass1.v + k2_v1 / 2, mass2.v + k2_v2 / 2);
    let k3_v2 = dt * accel2(mass1.the + k2_the1 / 2, mass2.the + k2_the2 / 2, mass1.v + k2_v1 / 2, mass2.v + k2_v2 / 2);
    let k3_the1 = dt * (mass1.v + k2_v1 / 2);
    let k3_the2 = dt * (mass2.v + k2_v2 / 2);

    let k4_v1 = dt * accel1(mass1.the + k3_the1, mass2.the + k3_the2, mass1.v + k3_v1, mass2.v + k3_v2);
    let k4_v2 = dt * accel2(mass1.the + k3_the1, mass2.the + k3_the2, mass1.v + k3_v1, mass2.v + k3_v2);
    let k4_the1 = dt * (mass1.v + k3_v1);
    let k4_the2 = dt * (mass2.v + k3_v2);

    mass1.the += (k1_the1 + 2 * k2_the1 + 2 * k3_the1 + k4_the1) / 6;
    mass2.the += (k1_the2 + 2 * k2_the2 + 2 * k3_the2 + k4_the2) / 6;
    mass1.v += (k1_v1 + 2 * k2_v1 + 2 * k3_v1 + k4_v1) / 6;
    mass2.v += (k1_v2 + 2 * k2_v2 + 2 * k3_v2 + k4_v2) / 6;

    coords1.x = pivot.x + length1 * Math.sin(mass1.the);
    coords1.y = pivot.y + length1 * Math.cos(mass1.the);
    coords2.x = coords1.x + length2 * Math.sin(mass2.the);
    coords2.y = coords1.y + length2 * Math.cos(mass2.the);
}

/** Renders the canvas every time step with the new coordinates */
function render() {
    // Sets the canvas background to black. If this wasn't here, it would create a trail effect
    // because the body background is also black
    ctx.fillStyle = 'black';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // Draws the two lines between the pivot and the pendulums
    drawLine(pivot, coords1);
    drawLine(coords1, coords2);

    // Draws the pivot and the pendulums
    drawCircle(pivot.x, pivot.y, 10, 'green');
    drawCircle(coords1.x, coords1.y, mass1.m / 6, 'red');
    drawCircle(coords2.x, coords2.y, mass2.m / 6, 'blue');
}

/** The main animation loop */
function animate() {
    setTimeout(function() {
        update(0.3); // The time step is 0.3 seconds
        render();
        requestAnimationFrame(animate);
    }, 1000 / 60);
}

requestAnimationFrame(animate);