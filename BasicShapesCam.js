//3456789_123456789_123456789_123456789_123456789_123456789_123456789_123456789_
// (JT: why the numbers? counts columns, helps me keep 80-char-wide listings)
//
// Chapter 5: ColoredTriangle.js (c) 2012 matsuda  AND
// Chapter 4: RotatingTriangle_withButtons.js (c) 2012 matsuda
// became:
//
// BasicShapes.js  MODIFIED for EECS 351-1, 
//									Northwestern Univ. Jack Tumblin
//		--converted from 2D to 4D (x,y,z,w) vertices
//		--extend to other attributes: color, surface normal, etc.
//		--demonstrate how to keep & use MULTIPLE colored shapes in just one
//			Vertex Buffer Object(VBO). 
//		--create several canonical 3D shapes borrowed from 'GLUT' library:
//		--Demonstrate how to make a 'stepped spiral' tri-strip,  and use it
//			to build a cylinder, sphere, and torus.
//
// Vertex shader program----------------------------------
var VSHADER_SOURCE =
	'uniform mat4 u_ModelMatrix;\n' +
	'uniform mat4 costum_Viewmatrix;\n' +
	'uniform mat4 costum_Projmatrix;\n' +
	'attribute vec4 a_Position;\n' +
	'attribute vec4 a_Color;\n' +
	'varying vec4 v_Color;\n' +
	'void main() {\n' +
	'  gl_Position = costum_Projmatrix * costum_Viewmatrix * u_ModelMatrix * a_Position;\n' +
	'  gl_PointSize = 10.0;\n' +
	'  v_Color = a_Color;\n' +
	'}\n';

// Fragment shader program----------------------------------
var FSHADER_SOURCE =
	//  '#ifdef GL_ES\n' +
	'precision mediump float;\n' +
	//  '#endif GL_ES\n' +
	'varying vec4 v_Color;\n' +
	'void main() {\n' +
	'  gl_FragColor = v_Color;\n' +
	'}\n';

// Global Variables
var ANGLE_STEP = 45.0;		// Rotation angle rate (degrees/second)
var floatsPerVertex = 7;	// # of Float32Array elements used for each vertex
// (x,y,z,w)position + (r,g,b)color
// Later, see if you can add:
// (x,y,z) surface normal + (tx,ty) texture addr.
var g_vertsMax = 0;                 // number of vertices held in the VBO 
// (global: replaces local 'n' variable)
var modelMatrix = new Matrix4(); 
var view = new Matrix4();
var projection = new Matrix4(); // Construct 4x4 matrix; contents get sent
// to the GPU/Shaders as a 'uniform' var.
var u_ModelMatrix;
var costum_Projmatrix;
var costum_Viewmatrix;  

//look_at variables
var x_location = 5
var y_location = 5
var z_location = 3

var x_refrence = 0
var y_refrence = 0
var z_refrence = 0

var x_up = 0
var y_up = 0
var z_up = 1


//------------For Animation---------------------------------------------
var g_isRun = true;                 // run/stop for animation; used in tick().
var g_lastMS = Date.now();    			// Timestamp for most-recently-drawn image; 
// in milliseconds; used by 'animate()' fcn 
// (now called 'timerAll()' ) to find time
// elapsed since last on-screen image.
var g_angle01 = 0;                  // initial rotation angle
var g_angle01Rate = 20.0;           // rotation speed, in degrees/second 

var g_angle02 = 0;                  // initial rotation angle
var g_angle02Rate = 15.0;

var g_angle03 = 0;
var g_angle03Rate = 200;

var g_angle04 = 0;
var g_angle04Rate = 200;

var g_angle05 = 0;
var g_angle05Rate = 50;

var w_is_pushed = false;
var s_is_pushed = false;

var w_last_pushed = true;

var ArrowDOWN_is_pushed = false;
var ArrowUP_is_pushed = false;
// rotation speed, in degrees/second 

//------------For mouse click-and-drag: -------------------------------
var g_isDrag = false;		// mouse-drag: true when user holds down mouse button
var g_xMclik = 0.0;			// last mouse button-down position (in CVV coords)
var g_yMclik = 0.0;
var g_xMdragTot = 0.0;	// total (accumulated) mouse-drag amounts (in CVV coords).
var g_yMdragTot = 0.0;
var g_digits = 5;			// DIAGNOSTICS: # of digits to print in console.log (
//    console.log('xVal:', xVal.toFixed(g_digits)); // print 5 digits

function main() {
	//==============================================================================
	// Retrieve <canvas> element
	var canvas = document.getElementById('webgl');

	// Get the rendering context for WebGL
	var gl = getWebGLContext(canvas);
	if (!gl) {
		console.log('Failed to get the rendering context for WebGL');
		return;
	}

	// Initialize shaders
	if (!initShaders(gl, VSHADER_SOURCE, FSHADER_SOURCE)) {
		console.log('Failed to intialize shaders.');
		return;
	}

	// 
	var n = initVertexBuffer(gl);
	if (n < 0) {
		console.log('Failed to set the vertex information');
		return;
	}

	// Specify the color for clearing <canvas>
	gl.clearColor(0.0, 0.0, 0.0, 1.0);

	// NEW!! Enable 3D depth-test when drawing: don't over-draw at any pixel 
	// unless the new Z value is closer to the eye than the old one..
	//	gl.depthFunc(gl.LESS);			 // WebGL default setting: (default)
	gl.enable(gl.DEPTH_TEST);

	window.addEventListener("keydown", myKeyDown, false);
	// After each 'keydown' event, call the 'myKeyDown()' function.  The 'false' 
	// arg (default) ensures myKeyDown() call in 'bubbling', not 'capture' stage)
	// ( https://www.w3schools.com/jsref/met_document_addeventlistener.asp )
	window.addEventListener("keyup", myKeyUp, false);
	// Called when user RELEASES the key.  Now rarely used...

	// MOUSE:
	// Create 'event listeners' for a few vital mouse events 
	// (others events are available too... google it!).  
	window.addEventListener("mousedown", myMouseDown);
	// (After each 'mousedown' event, browser calls the myMouseDown() fcn.)
	window.addEventListener("mousemove", myMouseMove);
	window.addEventListener("mouseup", myMouseUp);
	window.addEventListener("click", myMouseClick);
	window.addEventListener("dblclick", myMouseDblClick);

	//==============================================================================
	// STEP 4:   REMOVE This "reversed-depth correction"
	//       when you apply any of the 3D camera-lens transforms: 
	//      (e.g. Matrix4 member functions 'perspective(), frustum(), ortho() ...)
	//======================REVERSED-DEPTH Correction===============================



	//here changes
	 u_ModelMatrix = gl.getUniformLocation(gl.program, 'u_ModelMatrix');
	 costum_Viewmatrix = gl.getUniformLocation(gl.program, 'costum_Viewmatrix');
	 costum_Projmatrix = gl.getUniformLocation(gl.program, 'costum_Projmatrix');
	if (!u_ModelMatrix) {
		console.log('Failed to get the storage location of u_ModelMatrix');
		return;
	}

	if (!costum_Projmatrix || !costum_Viewmatrix) {

		console.log('failed to get stroage for costumes')
		return;
	}






	// Create, init current rotation angle value in JavaScript
	var currentAngle = 0.0;

	//-----------------  
	// Start drawing: create 'tick' variable whose value is this function:
	var tick = function () {
		currentAngle = animate(currentAngle);  // Update the rotation angle
		draw_two(gl, n, currentAngle, modelMatrix, u_ModelMatrix);   // Draw shapes
		// report current angle on console
		//console.log('currentAngle=',currentAngle);
		requestAnimationFrame(tick, canvas);
		// Request that the browser re-draw the webpage
	};
	tick();							// start (and continue) animation: draw current image

}

function initVertexBuffer(gl) {
	//==============================================================================
	// Create one giant vertex buffer object (VBO) that holds all vertices for all
	// shapes.

	// Make each 3D shape in its own array of vertices:
	make_old_shapes();				// create, fill the torVerts array
	makeGroundGrid();
	makeCylinder();
	makeTorus();
	make_Axes();				// create, fill the gndVerts array
	// how many floats total needed to store all shapes?
	var mySiz = ( gndVerts.length + old_shapes_vers.length + cylVerts.length + torVerts.length + axes_vars.length);

	// How many vertices total?
	var nn = mySiz / floatsPerVertex;
	console.log('nn is', nn, 'mySiz is', mySiz, 'floatsPerVertex is', floatsPerVertex);
	// Copy all shapes into one big Float32 array:
	var colorShapes = new Float32Array(mySiz);
	// Copy them:  remember where to start for each shape:
	old_ShapesStart = 0;							// we stored the cylinder first.
	for (i = 0, j = 0; j < old_shapes_vers.length; i++, j++) {
		colorShapes[i] = old_shapes_vers[j];
	}
	gndStart = i;							// we stored the cylinder first.
	for (j = 0; j < gndVerts.length; i++, j++) {
		colorShapes[i] = gndVerts[j];
	}
	cylStart = i;							// we stored the cylinder first.
	for (j = 0; j < cylVerts.length; i++, j++) {
		colorShapes[i] = cylVerts[j];
	}
	torStart = i;						// next, we'll store the torus;
	for (j = 0; j < torVerts.length; i++, j++) {
		colorShapes[i] = torVerts[j];
	}
	axesStart = i;
	for(j = 0; j < axes_vars.length; i++, j++){
		colorShapes[i] = axes_vars[j];
	}


	

	// Create a buffer object on the graphics hardware:
	var shapeBufferHandle = gl.createBuffer();
	if (!shapeBufferHandle) {
		console.log('Failed to create the shape buffer object');
		return false;
	}

	// Bind the the buffer object to target:
	gl.bindBuffer(gl.ARRAY_BUFFER, shapeBufferHandle);
	// Transfer data from Javascript array colorShapes to Graphics system VBO
	// (Use sparingly--may be slow if you transfer large shapes stored in files)
	gl.bufferData(gl.ARRAY_BUFFER, colorShapes, gl.STATIC_DRAW);

	//Get graphics system's handle for our Vertex Shader's position-input variable: 
	var a_Position = gl.getAttribLocation(gl.program, 'a_Position');
	if (a_Position < 0) {
		console.log('Failed to get the storage location of a_Position');
		return -1;
	}

	var FSIZE = colorShapes.BYTES_PER_ELEMENT; // how many bytes per stored value?

	// Use handle to specify how to retrieve **POSITION** data from our VBO:
	gl.vertexAttribPointer(
		a_Position, 	// choose Vertex Shader attribute to fill with data
		4, 						// how many values? 1,2,3 or 4.  (we're using x,y,z,w)
		gl.FLOAT, 		// data type for each value: usually gl.FLOAT
		false, 				// did we supply fixed-point data AND it needs normalizing?
		FSIZE * floatsPerVertex, // Stride -- how many bytes used to store each vertex?
		// (x,y,z,w, r,g,b) * bytes/value
		0);						// Offset -- now many bytes from START of buffer to the
	// value we will actually use?
	gl.enableVertexAttribArray(a_Position);
	// Enable assignment of vertex buffer object's position data

	// Get graphics system's handle for our Vertex Shader's color-input variable;
	var a_Color = gl.getAttribLocation(gl.program, 'a_Color');
	if (a_Color < 0) {
		console.log('Failed to get the storage location of a_Color');
		return -1;
	}
	// Use handle to specify how to retrieve **COLOR** data from our VBO:
	gl.vertexAttribPointer(
		a_Color, 				// choose Vertex Shader attribute to fill with data
		3, 							// how many values? 1,2,3 or 4. (we're using R,G,B)
		gl.FLOAT, 			// data type for each value: usually gl.FLOAT
		false, 					// did we supply fixed-point data AND it needs normalizing?
		FSIZE * 7, 			// Stride -- how many bytes used to store each vertex?
		// (x,y,z,w, r,g,b) * bytes/value
		FSIZE * 4);			// Offset -- how many bytes from START of buffer to the
	// value we will actually use?  Need to skip over x,y,z,w

	gl.enableVertexAttribArray(a_Color);
	// Enable assignment of vertex buffer object's position data

	//--------------------------------DONE!
	// Unbind the buffer object 
	gl.bindBuffer(gl.ARRAY_BUFFER, null);

	return nn;
}


function makeGroundGrid() {
	//==============================================================================
	// Create a list of vertices that create a large grid of lines in the x,y plane
	// centered at x=y=z=0.  Draw this shape using the GL_LINES primitive.

	var xcount = 100;			// # of lines to draw in x,y to make the grid.
	var ycount = 100;
	var xymax = 50.0;			// grid size; extends to cover +/-xymax in x and y.
	var xColr = new Float32Array([1.0, 1.0, 0.3]);	// bright yellow
	var yColr = new Float32Array([0.5, 1.0, 0.5]);	// bright green.

	// Create an (global) array to hold this ground-plane's vertices:
	gndVerts = new Float32Array(floatsPerVertex * 2 * (xcount + ycount));
	// draw a grid made of xcount+ycount lines; 2 vertices per line.

	var xgap = xymax / (xcount - 1);		// HALF-spacing between lines in x,y;
	var ygap = xymax / (ycount - 1);		// (why half? because v==(0line number/2))

	// First, step thru x values as we make vertical lines of constant-x:
	for (v = 0, j = 0; v < 2 * xcount; v++, j += floatsPerVertex) {
		if (v % 2 == 0) {	// put even-numbered vertices at (xnow, -xymax, 0)
			gndVerts[j] = -xymax + (v) * xgap;	// x
			gndVerts[j + 1] = -xymax;								// y
			gndVerts[j + 2] = 0.0;									// z
			gndVerts[j + 3] = 1.0;									// w.
		}
		else {				// put odd-numbered vertices at (xnow, +xymax, 0).
			gndVerts[j] = -xymax + (v - 1) * xgap;	// x
			gndVerts[j + 1] = xymax;								// y
			gndVerts[j + 2] = 0.0;									// z
			gndVerts[j + 3] = 1.0;									// w.
		}
		gndVerts[j + 4] = xColr[0];			// red
		gndVerts[j + 5] = xColr[1];			// grn
		gndVerts[j + 6] = xColr[2];			// blu
	}
	// Second, step thru y values as wqe make horizontal lines of constant-y:
	// (don't re-initialize j--we're adding more vertices to the array)
	for (v = 0; v < 2 * ycount; v++, j += floatsPerVertex) {
		if (v % 2 == 0) {		// put even-numbered vertices at (-xymax, ynow, 0)
			gndVerts[j] = -xymax;								// x
			gndVerts[j + 1] = -xymax + (v) * ygap;	// y
			gndVerts[j + 2] = 0.0;									// z
			gndVerts[j + 3] = 1.0;									// w.
		}
		else {					// put odd-numbered vertices at (+xymax, ynow, 0).
			gndVerts[j] = xymax;								// x
			gndVerts[j + 1] = -xymax + (v - 1) * ygap;	// y
			gndVerts[j + 2] = 0.0;									// z
			gndVerts[j + 3] = 1.0;									// w.
		}
		gndVerts[j + 4] = yColr[0];			// red
		gndVerts[j + 5] = yColr[1];			// grn
		gndVerts[j + 6] = yColr[2];			// blu
	}
}
function makeTorus() {
	//==============================================================================
	// 		Create a torus centered at the origin that circles the z axis.  
	// Terminology: imagine a torus as a flexible, cylinder-shaped bar or rod bent 
	// into a circle around the z-axis. The bent bar's centerline forms a circle
	// entirely in the z=0 plane, centered at the origin, with radius 'rbend'.  The 
	// bent-bar circle begins at (rbend,0,0), increases in +y direction to circle  
	// around the z-axis in counter-clockwise (CCW) direction, consistent with our
	// right-handed coordinate system.
	// 		This bent bar forms a torus because the bar itself has a circular cross-
	// section with radius 'rbar' and angle 'phi'. We measure phi in CCW direction 
	// around the bar's centerline, circling right-handed along the direction 
	// forward from the bar's start at theta=0 towards its end at theta=2PI.
	// 		THUS theta=0, phi=0 selects the torus surface point (rbend+rbar,0,0);
	// a slight increase in phi moves that point in -z direction and a slight
	// increase in theta moves that point in the +y direction.  
	// To construct the torus, begin with the circle at the start of the bar:
	//					xc = rbend + rbar*cos(phi); 
	//					yc = 0; 
	//					zc = -rbar*sin(phi);			(note negative sin(); right-handed phi)
	// and then rotate this circle around the z-axis by angle theta:
	//					x = xc*cos(theta) - yc*sin(theta) 	
	//					y = xc*sin(theta) + yc*cos(theta)
	//					z = zc
	// Simplify: yc==0, so
	//					x = (rbend + rbar*cos(phi))*cos(theta)
	//					y = (rbend + rbar*cos(phi))*sin(theta) 
	//					z = -rbar*sin(phi)
	// To construct a torus from a single triangle-strip, make a 'stepped spiral' 
	// along the length of the bent bar; successive rings of constant-theta, using 
	// the same design used for cylinder walls in 'makeCyl()' and for 'slices' in 
	// makeSphere().  Unlike the cylinder and sphere, we have no 'special case' 
	// for the first and last of these bar-encircling rings.
	//
	var rbend = 1.0;										// Radius of circle formed by torus' bent bar
	var rbar = 0.5;											// radius of the bar we bent to form torus
	var barSlices = 23;									// # of bar-segments in the torus: >=3 req'd;
	// more segments for more-circular torus
	var barSides = 13;										// # of sides of the bar (and thus the 
	// number of vertices in its cross-section)
	// >=3 req'd;
	// more sides for more-circular cross-section
	// for nice-looking torus with approx square facets, 
	//			--choose odd or prime#  for barSides, and
	//			--choose pdd or prime# for barSlices of approx. barSides *(rbend/rbar)
	// EXAMPLE: rbend = 1, rbar = 0.5, barSlices =23, barSides = 11.

	// Create a (global) array to hold this torus's vertices:
	torVerts = new Float32Array(floatsPerVertex * (2 * barSides * barSlices + 2));
	//	Each slice requires 2*barSides vertices, but 1st slice will skip its first 
	// triangle and last slice will skip its last triangle. To 'close' the torus,
	// repeat the first 2 vertices at the end of the triangle-strip.  Assume 7

	var phi = 0, theta = 0;										// begin torus at angles 0,0
	var thetaStep = 2 * Math.PI / barSlices;	// theta angle between each bar segment
	var phiHalfStep = Math.PI / barSides;		// half-phi angle between each side of bar
	// (WHY HALF? 2 vertices per step in phi)
	// s counts slices of the bar; v counts vertices within one slice; j counts
	// array elements (Float32) (vertices*#attribs/vertex) put in torVerts array.
	for (s = 0, j = 0; s < barSlices; s++) {		// for each 'slice' or 'ring' of the torus:
		for (v = 0; v < 2 * barSides; v++, j += 7) {		// for each vertex in this slice:
			if (v % 2 == 0) {	// even #'d vertices at bottom of slice,
				torVerts[j] = (rbend + rbar * Math.cos((v) * phiHalfStep)) *
					Math.cos((s) * thetaStep);
				//	x = (rbend + rbar*cos(phi)) * cos(theta)
				torVerts[j + 1] = (rbend + rbar * Math.cos((v) * phiHalfStep)) *
					Math.sin((s) * thetaStep);
				//  y = (rbend + rbar*cos(phi)) * sin(theta) 
				torVerts[j + 2] = -rbar * Math.sin((v) * phiHalfStep);
				//  z = -rbar  *   sin(phi)
				torVerts[j + 3] = 1.0;		// w
			}
			else {				// odd #'d vertices at top of slice (s+1);
				// at same phi used at bottom of slice (v-1)
				torVerts[j] = (rbend + rbar * Math.cos((v - 1) * phiHalfStep)) *
					Math.cos((s + 1) * thetaStep);
				//	x = (rbend + rbar*cos(phi)) * cos(theta)
				torVerts[j + 1] = (rbend + rbar * Math.cos((v - 1) * phiHalfStep)) *
					Math.sin((s + 1) * thetaStep);
				//  y = (rbend + rbar*cos(phi)) * sin(theta) 
				torVerts[j + 2] = -rbar * Math.sin((v - 1) * phiHalfStep);
				//  z = -rbar  *   sin(phi)
				torVerts[j + 3] = 1.0;		// w
			}
			torVerts[j + 4] = Math.random();		// random color 0.0 <= R < 1.0
			torVerts[j + 5] = Math.random();		// random color 0.0 <= G < 1.0
			torVerts[j + 6] = Math.random();		// random color 0.0 <= B < 1.0
		}
	}
	// Repeat the 1st 2 vertices of the triangle strip to complete the torus:
	torVerts[j] = rbend + rbar;	// copy vertex zero;
	//	x = (rbend + rbar*cos(phi==0)) * cos(theta==0)
	torVerts[j + 1] = 0.0;
	//  y = (rbend + rbar*cos(phi==0)) * sin(theta==0) 
	torVerts[j + 2] = 0.0;
	//  z = -rbar  *   sin(phi==0)
	torVerts[j + 3] = 1.0;		// w
	torVerts[j + 4] = Math.random();		// random color 0.0 <= R < 1.0
	torVerts[j + 5] = Math.random();		// random color 0.0 <= G < 1.0
	torVerts[j + 6] = Math.random();		// random color 0.0 <= B < 1.0
	j += 7; // go to next vertex:
	torVerts[j] = (rbend + rbar) * Math.cos(thetaStep);
	//	x = (rbend + rbar*cos(phi==0)) * cos(theta==thetaStep)
	torVerts[j + 1] = (rbend + rbar) * Math.sin(thetaStep);
	//  y = (rbend + rbar*cos(phi==0)) * sin(theta==thetaStep) 
	torVerts[j + 2] = 0.0;
	//  z = -rbar  *   sin(phi==0)
	torVerts[j + 3] = 1.0;		// w
	torVerts[j + 4] = Math.random();		// random color 0.0 <= R < 1.0
	torVerts[j + 5] = Math.random();		// random color 0.0 <= G < 1.0
	torVerts[j + 6] = Math.random();		// random color 0.0 <= B < 1.0
}
function makeCylinder() {
	//==============================================================================
	// Make a cylinder shape from one TRIANGLE_STRIP drawing primitive, using the
	// 'stepped spiral' design described in notes.
	// Cylinder center at origin, encircles z axis, radius 1, top/bottom at z= +/-1.
	//
	var ctrColr = new Float32Array([0.2, 0.2, 0.2]);	// dark gray
	var topColr = new Float32Array([0.4, 0.7, 0.4]);	// light green
	var botColr = new Float32Array([0.5, 0.5, 1.0]);	// light blue
	var capVerts = 16;	// # of vertices around the topmost 'cap' of the shape
	var botRadius = 1.6;		// radius of bottom of cylinder (top always 1.0)

	// Create a (global) array to hold this cylinder's vertices;
	cylVerts = new Float32Array(((capVerts * 6) - 2) * floatsPerVertex);
	// # of vertices * # of elements needed to store them. 

	// Create circle-shaped top cap of cylinder at z=+1.0, radius 1.0
	// v counts vertices: j counts array elements (vertices * elements per vertex)
	for (v = 1, j = 0; v < 2 * capVerts; v++, j += floatsPerVertex) {
		// skip the first vertex--not needed.
		if (v % 2 == 0) {				// put even# vertices at center of cylinder's top cap:
			cylVerts[j] = 0.0; 			// x,y,z,w == 0,0,1,1
			cylVerts[j + 1] = 0.0;
			cylVerts[j + 2] = 1.0;
			cylVerts[j + 3] = 1.0;			// r,g,b = topColr[]
			cylVerts[j + 4] = ctrColr[0];
			cylVerts[j + 5] = ctrColr[1];
			cylVerts[j + 6] = ctrColr[2];
		}
		else { 	// put odd# vertices around the top cap's outer edge;
			// x,y,z,w == cos(theta),sin(theta), 1.0, 1.0
			// 					theta = 2*PI*((v-1)/2)/capVerts = PI*(v-1)/capVerts
			cylVerts[j] = Math.cos(Math.PI * (v - 1) / capVerts);			// x
			cylVerts[j + 1] = Math.sin(Math.PI * (v - 1) / capVerts);			// y
			//	(Why not 2*PI? because 0 < =v < 2*capVerts, so we
			//	 can simplify cos(2*PI * (v-1)/(2*capVerts))
			cylVerts[j + 2] = 1.0;	// z
			cylVerts[j + 3] = 1.0;	// w.
			// r,g,b = topColr[]
			cylVerts[j + 4] = topColr[0];
			cylVerts[j + 5] = topColr[1];
			cylVerts[j + 6] = topColr[2];
		}
	}
	// Create the cylinder side walls, made of 2*capVerts vertices.
	// v counts vertices within the wall; j continues to count array elements
	for (v = 0; v < 2 * capVerts; v++, j += floatsPerVertex) {
		if (v % 2 == 0)	// position all even# vertices along top cap:
		{
			cylVerts[j] = Math.cos(Math.PI * (v) / capVerts);		// x
			cylVerts[j + 1] = Math.sin(Math.PI * (v) / capVerts);		// y
			cylVerts[j + 2] = 1.0;	// z
			cylVerts[j + 3] = 1.0;	// w.
			// r,g,b = topColr[]
			cylVerts[j + 4] = topColr[0];
			cylVerts[j + 5] = topColr[1];
			cylVerts[j + 6] = topColr[2];
		}
		else		// position all odd# vertices along the bottom cap:
		{
			cylVerts[j] = botRadius * Math.cos(Math.PI * (v - 1) / capVerts);		// x
			cylVerts[j + 1] = botRadius * Math.sin(Math.PI * (v - 1) / capVerts);		// y
			cylVerts[j + 2] = -1.0;	// z
			cylVerts[j + 3] = 1.0;	// w.
			// r,g,b = topColr[]
			cylVerts[j + 4] = botColr[0];
			cylVerts[j + 5] = botColr[1];
			cylVerts[j + 6] = botColr[2];
		}
	}
	// Create the cylinder bottom cap, made of 2*capVerts -1 vertices.
	// v counts the vertices in the cap; j continues to count array elements
	for (v = 0; v < (2 * capVerts - 1); v++, j += floatsPerVertex) {
		if (v % 2 == 0) {	// position even #'d vertices around bot cap's outer edge
			cylVerts[j] = botRadius * Math.cos(Math.PI * (v) / capVerts);		// x
			cylVerts[j + 1] = botRadius * Math.sin(Math.PI * (v) / capVerts);		// y
			cylVerts[j + 2] = -1.0;	// z
			cylVerts[j + 3] = 1.0;	// w.
			// r,g,b = topColr[]
			cylVerts[j + 4] = botColr[0];
			cylVerts[j + 5] = botColr[1];
			cylVerts[j + 6] = botColr[2];
		}
		else {				// position odd#'d vertices at center of the bottom cap:
			cylVerts[j] = 0.0; 			// x,y,z,w == 0,0,-1,1
			cylVerts[j + 1] = 0.0;
			cylVerts[j + 2] = -1.0;
			cylVerts[j + 3] = 1.0;			// r,g,b = botColr[]
			cylVerts[j + 4] = botColr[0];
			cylVerts[j + 5] = botColr[1];
			cylVerts[j + 6] = botColr[2];
		}
	}
}

function make_Axes(){

 axes_vars = new Float32Array([

0,0,0,1,1,0,0,
1,0,0,1,1,0,0,

0,0,0,1,0,1,0,
0,1,0,1,0,1,0,

0,0,0,1,0,0,1,
0,0,1,1,0,0,1




 ])





}

//add old shapes

function make_old_shapes(){

	//head variables
	//face
	var x1 = 0; var y1 = 0; var z1 = 0;
	var x2 = 0.3; var y2 = 0.5; var z2 = 0.2;

	var c1 = 1; var c2 = 1; var c3 = 0.7;
	// cap variables

	//lower cap
	var x11 = -0.15; var y11 = 0.5; var z11 = -0.15;
	var x12 = 0.45; var y12 = 0.6; var z12 = 0.35;

	var c11 = 0; var c12 = 0; var c13 = 0;

	//uper cap
	var x21 = 0; var y21 = 0.6; var z21 = 0;
	var x22 = 0.3; var y22 = 1; var z22 = 0.2;

	var c21 = 0; var c22 = 0; var c23 = 0;

	//head end

	//arm variables
	//arm

	var x31 = 0; var y31 = 0; var z31 = 0;
	var x32 = 0.2; var y32 = 0.5; var z32 = 0.2;

	var c31 = 1; var c32 = 1; var c33 = 0.7;

	//hand
	var x41 = -0.025; var y41 = -0.3; var z41 = -0.025;
	var x42 = 0.225; var y42 = 0; var z42 = 0.225;

	var c41 = 0.647059; var c42 = 0.164706; var c43 = 0.164706;

	//body variables
	//body

	var x51 = 0; var y51 = 0; var z51 = 0;
	var x52 = 0.4; var y52 = 0.8; var z52 = 0.2;

	var c51 = 0; var c52 = 1; var c53 = 0.7;

	//leg variables
	//upper part

	var x61 = 0; var y61 = 0; var z61 = 0;
	var x62 = 0.2; var y62 = 0.7; var z62 = 0.2;

	var c61 = 1; var c62 = 1; var c63 = 0.7;

	//foot
	var x71 = -0.0; var y71 = -0.2; var z71 = -0.3;
	var x72 = 0.2; var y72 = 0; var z72 = 0.225;

	var c71 = 0; var c72 = 0; var c73 = 0;

	//cowboy end

	//dog start
	// head variables

	//face
	var x81 = 0; var y81 = 0; var z81 = 0;
	var x82 = 0.3; var y82 = 0.3; var z82 = 0.2;

	var c81 = 1; var c82 = 0.8; var c83 = 0.7;

	//left ear
	var x91 = 0; var y91 = 0.3; var z91 = 0;
	var x92 = 0.1; var y92 = 0.4; var z92 = 0.1;

	var c91 = 0.8; var c92 = 0.8; var c93 = 0.7;
	//right ear
	var xa1 = 0.2; var ya1 = 0.3; var za1 = 0;
	var xa2 = 0.3; var ya2 = 0.4; var za2 = 0.1;

	var ca1 = 0.8; var ca2 = 0.8; var ca3 = 0.7;
	//head end
	//body start

	var xb1 = 0; var yb1 = 0; var zb1 = 0;
	var xb2 = 0.4; var yb2 = 0.2; var zb2 = 0.6;

	var cb1 = 0.9; var cb2 = 0.8; var cb3 = 0.7;
	//body end

	//legs start
	var xc1 = 0; var yc1 = 0; var zc1 = 0;
	var xc2 = 0.2; var yc2 = 0.35; var zc2 = 0.2;

	var cc1 = 1; var cc2 = 0.8; var cc3 = 0.7;

	//paw
	var xd1 = -0.025; var yd1 = -0.2; var zd1 = -0.025;
	var xd2 = 0.225; var yd2 = 0; var zd2 = 0.225;

	var cd1 = 0.8; var cd2 = 0.8; var cd3 = 0.7;

	//legs end
	//pig end

	old_shapes_vers =  new Float32Array([
		// Vertex coordinates(x,y,z,w) and color (R,G,B) for a color tetrahedron:
		//		Apex on +z axis; equilateral triangle base at z=0
		/*	Nodes:
				 0.0,	 0.0, sq2, 1.0,			1.0, 	1.0,	1.0,	// Node 0 (apex, +z axis;  white)
			 c30, -0.5, 0.0, 1.0, 		0.0,  0.0,  1.0, 	// Node 1 (base: lower rt; red)
			 0.0,  1.0, 0.0, 1.0,  		1.0,  0.0,  0.0,	// Node 2 (base: +y axis;  grn)
			-c30, -0.5, 0.0, 1.0, 		0.0,  1.0,  0.0, 	// Node 3 (base:lower lft; blue)
		
		*/





		// HEAAAD PART START

		//face
		// Face : back
		x1, y1, z1, 1.0, c1, c2, c3,
		x2, y1, z1, 1.0, c1, c2, c3,
		x1, y2, z1, 1.0, c1, c2, c3,

		x2, y2, z1, 1.00, c1, c2, c3,
		x2, y1, z1, 1.00, c1, c2, c3,
		x1, y2, z1, 1.00, c1, c2, c3,

		// Face : front
		x1, y1, z2, 1.0, c1, c2, c3,
		x2, y1, z2, 1.0, c1, c2, c3,
		x1, y2, z2, 1.0, c1, c2, c3,

		x2, y2, z2, 1.00, c1, c2, c3,
		x2, y1, z2, 1.00, c1, c2, c3,
		x1, y2, z2, 1.00, c1, c2, c3,

		//right

		x2, y1, z2, 1.0, c1, c2, c3,
		x2, y1, z1, 1.0, c1, c2, c3,
		x2, y2, z2, 1.00, c1, c2, c3,


		x2, y2, z1, 1.00, c1, c2, c3,
		x2, y2, z2, 1.00, c1, c2, c3,
		x2, y1, z1, 1.00, c1, c2, c3,

		//left


		x1, y1, z2, 1.0, c1, c2, c3,
		x1, y1, z1, 1.0, c1, c2, c3,
		x1, y2, z2, 1.00, c1, c2, c3,


		x1, y2, z1, 1.00, c1, c2, c3,
		x1, y2, z2, 1.00, c1, c2, c3,
		x1, y1, z1, 1.00, c1, c2, c3,

		//bottom
		x1, y1, z1, 1.0, c1, c2, c3,
		x1, y1, z2, 1.0, c1, c2, c3,
		x2, y1, z1, 1.0, c1, c2, c3,

		x2, y1, z2, 1.0, c1, c2, c3,
		x2, y1, z1, 1.0, c1, c2, c3,
		x1, y1, z2, 1.0, c1, c2, c3,

		//top

		x1, y2, z1, 1.0, c1, c2, c3,
		x1, y2, z2, 1.0, c1, c2, c3,
		x2, y2, z1, 1.0, c1, c2, c3,

		x2, y2, z2, 1.0, c1, c2, c3,
		x2, y2, z1, 1.0, c1, c2, c3,
		x1, y2, z2, 1.0, c1, c2, c3,

		// cap lower part

		//back
		x11, y11, z11, 1.0, c11, c12, c13,
		x12, y11, z11, 1.0, c11, c12, c13,
		x11, y12, z11, 1.0, c11, c12, c13,

		x12, y12, z11, 1.00, c11, c12, c13,
		x12, y11, z11, 1.00, c11, c12, c13,
		x11, y12, z11, 1.00, c11, c12, c13,

		// Face : front
		x11, y11, z12, 1.0, c11, c12, c13,
		x12, y11, z12, 1.0, c11, c12, c13,
		x11, y12, z12, 1.0, c11, c12, c13,

		x12, y12, z12, 1.00, c11, c12, c13,
		x12, y11, z12, 1.00, c11, c12, c13,
		x11, y12, z12, 1.00, c11, c12, c13,

		//right

		x12, y11, z12, 1.0, c11, c12, c13,
		x12, y11, z11, 1.0, c11, c12, c13,
		x12, y12, z12, 1.00, c11, c12, c13,


		x12, y12, z11, 1.00, c11, c12, c13,
		x12, y12, z12, 1.00, c11, c12, c13,
		x12, y11, z11, 1.00, c11, c12, c13,

		//left


		x11, y11, z12, 1.0, c11, c12, c13,
		x11, y11, z11, 1.0, c11, c12, c13,
		x11, y12, z12, 1.00, c11, c12, c13,


		x11, y12, z11, 1.00, c11, c12, c13,
		x11, y12, z12, 1.00, c11, c12, c13,
		x11, y11, z11, 1.00, c11, c12, c13,

		//bottom
		x11, y11, z11, 1.0, c11, c12, c13,
		x11, y11, z12, 1.0, c11, c12, c13,
		x12, y11, z11, 1.0, c11, c12, c13,

		x12, y11, z12, 1.0, c11, c12, c13,
		x12, y11, z11, 1.0, c11, c12, c13,
		x11, y11, z12, 1.0, c11, c12, c13,

		//top

		x11, y12, z11, 1.0, c11, c12, c13,
		x11, y12, z12, 1.0, c11, c12, c13,
		x12, y12, z11, 1.0, c11, c12, c13,

		x12, y12, z12, 1.0, c11, c12, c13,
		x12, y12, z11, 1.0, c11, c12, c13,
		x11, y12, z12, 1.0, c11, c12, c13,



		//upper part

		//back
		x21, y21, z21, 1.0, c21, c22, c23,
		x22, y21, z21, 1.0, c21, c22, c23,
		x21, y22, z21, 1.0, c21, c22, c23,
		x22, y22, z21, 1.0, c21, c22, c23,
		x22, y21, z21, 1.0, c21, c22, c23,
		x21, y22, z21, 1.0, c21, c22, c23,

		// Face : front
		x21, y21, z22, 1.0, c21, c22, c23,
		x22, y21, z22, 1.0, c21, c22, c23,
		x21, y22, z22, 1.0, c21, c22, c23,
		x22, y22, z22, 1.0, c21, c22, c23,
		x22, y21, z22, 1.0, c21, c22, c23,
		x21, y22, z22, 1.0, c21, c22, c23,

		//right

		x22, y21, z22, 1.0, c21, c22, c23,
		x22, y21, z21, 1.0, c21, c22, c23,
		x22, y22, z22, 1.0, c21, c22, c23,
		x22, y22, z21, 1.0, c21, c22, c23,
		x22, y22, z22, 1.0, c21, c22, c23,
		x22, y21, z21, 1.0, c21, c22, c23,

		//left


		x21, y21, z22, 1.0, c21, c22, c23,
		x21, y21, z21, 1.0, c21, c22, c23,
		x21, y22, z22, 1.0, c21, c22, c23,
		x21, y22, z21, 1.0, c21, c22, c23,
		x21, y22, z22, 1.0, c21, c22, c23,
		x21, y21, z21, 1.0, c21, c22, c23,

		//bottom
		x21, y21, z21, 1.0, c21, c22, c23,
		x21, y21, z22, 1.0, c21, c22, c23,
		x22, y21, z21, 1.0, c21, c22, c23,
		x22, y21, z22, 1.0, c21, c22, c23,
		x22, y21, z21, 1.0, c21, c22, c23,
		x21, y21, z22, 1.0, c21, c22, c23,

		//top

		x21, y22, z21, 1.0, 1, 0, 0,
		x21, y22, z22, 1.0, 0, 1, 0,
		x22, y22, z21, 1.0, 0, 0, 1,
		x22, y22, z22, 1.0, c21, c22, c23,
		x22, y22, z21, 1.0, c21, c22, c23,
		x21, y22, z22, 1.0, c21, c22, c23,

		// HEADPART END

		//ARMPART START
		//arm
		//arm : back
		x31, y31, z31, 1.0, c31, c32, c33,
		x32, y31, z31, 1.0, c31, c32, c33,
		x31, y32, z31, 1.0, c31, c32, c33,
		x32, y32, z31, 1.0, c31, c32, c33,
		x32, y31, z31, 1.0, c31, c32, c33,
		x31, y32, z31, 1.0, c31, c32, c33,

		//Face : front
		x31, y31, z32, 1.0, c31, c32, c33,
		x32, y31, z32, 1.0, c31, c32, c33,
		x31, y32, z32, 1.0, c31, c32, c33,
		x32, y32, z32, 1.0, c31, c32, c33,
		x32, y31, z32, 1.0, c31, c32, c33,
		x31, y32, z32, 1.0, c31, c32, c33,

		//right

		x32, y31, z32, 1.0, c31, c32, c33,
		x32, y31, z31, 1.0, c31, c32, c33,
		x32, y32, z32, 1.0, c31, c32, c33,
		x32, y32, z31, 1.0, c31, c32, c33,
		x32, y32, z32, 1.0, c31, c32, c33,
		x32, y31, z31, 1.0, c31, c32, c33,

		//left


		x31, y31, z32, 1.0, c31, c32, c33,
		x31, y31, z31, 1.0, c31, c32, c33,
		x31, y32, z32, 1.0, c31, c32, c33,
		x31, y32, z31, 1.0, c31, c32, c33,
		x31, y32, z32, 1.0, c31, c32, c33,
		x31, y31, z31, 1.0, c31, c32, c33,

		//	bottom
		x31, y31, z31, 1.0, c31, c32, c33,
		x31, y31, z32, 1.0, c31, c32, c33,
		x32, y31, z31, 1.0, c31, c32, c33,
		x32, y31, z32, 1.0, c31, c32, c33,
		x32, y31, z31, 1.0, c31, c32, c33,
		x31, y31, z32, 1.0, c31, c32, c33,

		//	top

		x31, y32, z31, 1.0, c31, c32, c33,
		x31, y32, z32, 1.0, c31, c32, c33,
		x32, y32, z31, 1.0, c31, c32, c33,
		x32, y32, z32, 1.0, c31, c32, c33,
		x32, y32, z31, 1.0, c31, c32, c33,
		x31, y32, z32, 1.0, c31, c32, c33,

		//Hand part
		//back
		x41, y41, z41, 1.0, c41, c42, c43,
		x42, y41, z41, 1.0, c41, c42, c43,
		x41, y42, z41, 1.0, c41, c42, c43,
		x42, y42, z41, 1.0, c41, c42, c43,
		x42, y41, z41, 1.0, c41, c42, c43,
		x41, y42, z41, 1.0, c41, c42, c43,

		//Face : front
		x41, y41, z42, 1.0, c41, c42, c43,
		x42, y41, z42, 1.0, c41, c42, c43,
		x41, y42, z42, 1.0, c41, c42, c43,
		x42, y42, z42, 1.0, c41, c42, c43,
		x42, y41, z42, 1.0, c41, c42, c43,
		x41, y42, z42, 1.0, c41, c42, c43,

		//	right

		x42, y41, z42, 1.0, c41, c42, c43,
		x42, y41, z41, 1.0, c41, c42, c43,
		x42, y42, z42, 1.0, c41, c42, c43,
		x42, y42, z41, 1.0, c41, c42, c43,
		x42, y42, z42, 1.0, c41, c42, c43,
		x42, y41, z41, 1.0, c41, c42, c43,

		//left


		x41, y41, z42, 1.0, c41, c42, c43,
		x41, y41, z41, 1.0, c41, c42, c43,
		x41, y42, z42, 1.0, c41, c42, c43,
		x41, y42, z41, 1.0, c41, c42, c43,
		x41, y42, z42, 1.0, c41, c42, c43,
		x41, y41, z41, 1.0, c41, c42, c43,

		//	bottom
		x41, y41, z41, 1.0, c41, c42, c43,
		x41, y41, z42, 1.0, c41, c42, c43,
		x42, y41, z41, 1.0, c41, c42, c43,
		x42, y41, z42, 1.0, c41, c42, c43,
		x42, y41, z41, 1.0, c41, c42, c43,
		x41, y41, z42, 1.0, c41, c42, c43,

		//	top

		x41, y42, z41, 1.0, c41, c42, c43,
		x41, y42, z42, 1.0, c41, c42, c43,
		x42, y42, z41, 1.0, c41, c42, c43,
		x42, y42, z42, 1.0, c41, c42, c43,
		x42, y42, z41, 1.0, c41, c42, c43,
		x41, y42, z42, 1.0, c41, c42, c43,

		//body start

		//back
		x51, y51, z51, 1.0, c51, c52, c53,
		x52, y51, z51, 1.0, c51, c52, c53,
		x51, y52, z51, 1.0, c51, c52, c53,
		x52, y52, z51, 1.0, c51, c52, c53,
		x52, y51, z51, 1.0, c51, c52, c53,
		x51, y52, z51, 1.0, c51, c52, c53,

		// Face : front
		x51, y51, z52, 1.0, c51, c52, c53,
		x52, y51, z52, 1.0, c51, c52, c53,
		x51, y52, z52, 1.0, c51, c52, c53,
		x52, y52, z52, 1.0, c51, c52, c53,
		x52, y51, z52, 1.0, c51, c52, c53,
		x51, y52, z52, 1.0, c51, c52, c53,

		//right

		x52, y51, z52, 1.0, c51, c52, c53,
		x52, y51, z51, 1.0, c51, c52, c53,
		x52, y52, z52, 1.0, c51, c52, c53,
		x52, y52, z51, 1.0, c51, c52, c53,
		x52, y52, z52, 1.0, c51, c52, c53,
		x52, y51, z51, 1.0, c51, c52, c53,

		//left


		x51, y51, z52, 1.0, c51, c52, c53,
		x51, y51, z51, 1.0, c51, c52, c53,
		x51, y52, z52, 1.0, c51, c52, c53,
		x51, y52, z51, 1.0, c51, c52, c53,
		x51, y52, z52, 1.0, c51, c52, c53,
		x51, y51, z51, 1.0, c51, c52, c53,

		//bottom
		x51, y51, z51, 1.0, c51, c52, c53,
		x51, y51, z52, 1.0, c51, c52, c53,
		x52, y51, z51, 1.0, c51, c52, c53,
		x52, y51, z52, 1.0, c51, c52, c53,
		x52, y51, z51, 1.0, c51, c52, c53,
		x51, y51, z52, 1.0, c51, c52, c53,

		//top

		x51, y52, z51, 1.0, c51, c52, c53,
		x51, y52, z52, 1.0, c51, c52, c53,
		x52, y52, z51, 1.0, c51, c52, c53,
		x52, y52, z52, 1.0, c51, c52, c53,
		x52, y52, z51, 1.0, c51, c52, c53,
		x51, y52, z52, 1.0, c51, c52, c53,

		// body end

		//leg start
		//leg uper part
		x61, y61, z61, 1.0, c61, c62, c63,
		x62, y61, z61, 1.0, c61, c62, c63,
		x61, y62, z61, 1.0, c61, c62, c63,
		x62, y62, z61, 1.0, c61, c62, c63,
		x62, y61, z61, 1.0, c61, c62, c63,
		x61, y62, z61, 1.0, c61, c62, c63,

		//Face : front
		x61, y61, z62, 1.0, c61, c62, c63,
		x62, y61, z62, 1.0, c61, c62, c63,
		x61, y62, z62, 1.0, c61, c62, c63,
		x62, y62, z62, 1.0, c61, c62, c63,
		x62, y61, z62, 1.0, c61, c62, c63,
		x61, y62, z62, 1.0, c61, c62, c63,

		//right

		x62, y61, z62, 1.0, c61, c62, c63,
		x62, y61, z61, 1.0, c61, c62, c63,
		x62, y62, z62, 1.0, c61, c62, c63,
		x62, y62, z61, 1.0, c61, c62, c63,
		x62, y62, z62, 1.0, c61, c62, c63,
		x62, y61, z61, 1.0, c61, c62, c63,

		//left


		x61, y61, z62, 1.0, c61, c62, c63,
		x61, y61, z61, 1.0, c61, c62, c63,
		x61, y62, z62, 1.0, c61, c62, c63,
		x61, y62, z61, 1.0, c61, c62, c63,
		x61, y62, z62, 1.0, c61, c62, c63,
		x61, y61, z61, 1.0, c61, c62, c63,

		//	bottom
		x61, y61, z61, 1.0, c61, c62, c63,
		x61, y61, z62, 1.0, c61, c62, c63,
		x62, y61, z61, 1.0, c61, c62, c63,
		x62, y61, z62, 1.0, c61, c62, c63,
		x62, y61, z61, 1.0, c61, c62, c63,
		x61, y61, z62, 1.0, c61, c62, c63,

		//	top

		x61, y62, z61, 1.0, c61, c62, c63,
		x61, y62, z62, 1.0, c61, c62, c63,
		x62, y62, z61, 1.0, c61, c62, c63,
		x62, y62, z62, 1.0, c61, c62, c63,
		x62, y62, z61, 1.0, c61, c62, c63,
		x61, y62, z62, 1.0, c61, c62, c63,

		//Hand part
		//back
		x71, y71, z71, 1.0, c71, c72, c73,
		x72, y71, z71, 1.0, c71, c72, c73,
		x71, y72, z71, 1.0, c71, c72, c73,
		x72, y72, z71, 1.0, c71, c72, c73,
		x72, y71, z71, 1.0, c71, c72, c73,
		x71, y72, z71, 1.0, c71, c72, c73,

		//Face : front
		x71, y71, z72, 1.0, c71, c72, c73,
		x72, y71, z72, 1.0, c71, c72, c73,
		x71, y72, z72, 1.0, c71, c72, c73,
		x72, y72, z72, 1.0, c71, c72, c73,
		x72, y71, z72, 1.0, c71, c72, c73,
		x71, y72, z72, 1.0, c71, c72, c73,

		//	right

		x72, y71, z72, 1.0, c71, c72, c73,
		x72, y71, z71, 1.0, c71, c72, c73,
		x72, y72, z72, 1.0, c71, c72, c73,
		x72, y72, z71, 1.0, c71, c72, c73,
		x72, y72, z72, 1.0, c71, c72, c73,
		x72, y71, z71, 1.0, c71, c72, c73,

		//left


		x71, y71, z72, 1.0, c71, c72, c73,
		x71, y71, z71, 1.0, c71, c72, c73,
		x71, y72, z72, 1.0, c71, c72, c73,
		x71, y72, z71, 1.0, c71, c72, c73,
		x71, y72, z72, 1.0, c71, c72, c73,
		x71, y71, z71, 1.0, c71, c72, c73,

		//	bottom
		x71, y71, z71, 1.0, c71, c72, c73,
		x71, y71, z72, 1.0, c71, c72, c73,
		x72, y71, z71, 1.0, c71, c72, c73,
		x72, y71, z72, 1.0, c71, c72, c73,
		x72, y71, z71, 1.0, c71, c72, c73,
		x71, y71, z72, 1.0, c71, c72, c73,

		//	top

		x71, y72, z71, 1.0, c71, c72, c73,
		x71, y72, z72, 1.0, c71, c72, c73,
		x72, y72, z71, 1.0, c71, c72, c73,
		x72, y72, z72, 1.0, c71, c72, c73,
		x72, y72, z71, 1.0, c71, c72, c73,
		x71, y72, z72, 1.0, c71, c72, c73,
		//leg end

		//end cowboy
		//	dog	start

		//face


		x81, y81, z81, 1.0, c81, c82, c83,
		x82, y81, z81, 1.0, c81, c82, c83,
		x81, y82, z81, 1.0, c81, c82, c83,
		x82, y82, z81, 1.0, c81, c82, c83,
		x82, y81, z81, 1.0, c81, c82, c83,
		x81, y82, z81, 1.0, c81, c82, c83,

		//Face : front
		x81, y81, z82, 1.0, c81, c82, c83,
		x82, y81, z82, 1.0, c81, c82, c83,
		x81, y82, z82, 1.0, c81, c82, c83,
		x82, y82, z82, 1.0, c81, c82, c83,
		x82, y81, z82, 1.0, c81, c82, c83,
		x81, y82, z82, 1.0, c81, c82, c83,

		//	right

		x82, y81, z82, 1.0, c81, c82, c83,
		x82, y81, z81, 1.0, c81, c82, c83,
		x82, y82, z82, 1.0, c81, c82, c83,
		x82, y82, z81, 1.0, c81, c82, c83,
		x82, y82, z82, 1.0, c81, c82, c83,
		x82, y81, z81, 1.0, c81, c82, c83,

		//left


		x81, y81, z82, 1.0, c81, c82, c83,
		x81, y81, z81, 1.0, c81, c82, c83,
		x81, y82, z82, 1.0, c81, c82, c83,
		x81, y82, z81, 1.0, c81, c82, c83,
		x81, y82, z82, 1.0, c81, c82, c83,
		x81, y81, z81, 1.0, c81, c82, c83,

		//	bottom
		x81, y81, z81, 1.0, c81, c82, c83,
		x81, y81, z82, 1.0, c81, c82, c83,
		x82, y81, z81, 1.0, c81, c82, c83,
		x82, y81, z82, 1.0, c81, c82, c83,
		x82, y81, z81, 1.0, c81, c82, c83,
		x81, y81, z82, 1.0, c81, c82, c83,

		//	top

		x81, y82, z81, 1.0, c81, c82, c83,
		x81, y82, z82, 1.0, c81, c82, c83,
		x82, y82, z81, 1.0, c81, c82, c83,
		x82, y82, z82, 1.0, c81, c82, c83,
		x82, y82, z81, 1.0, c81, c82, c83,
		x81, y82, z82, 1.0, c81, c82, c83,

		//left ear


		x91, y91, z91, 1.0, c91, c92, c93,
		x92, y91, z91, 1.0, c91, c92, c93,
		x91, y92, z91, 1.0, c91, c92, c93,
		x92, y92, z91, 1.0, c91, c92, c93,
		x92, y91, z91, 1.0, c91, c92, c93,
		x91, y92, z91, 1.0, c91, c92, c93,

		//Face : front
		x91, y91, z92, 1.0, c91, c92, c93,
		x92, y91, z92, 1.0, c91, c92, c93,
		x91, y92, z92, 1.0, c91, c92, c93,
		x92, y92, z92, 1.0, c91, c92, c93,
		x92, y91, z92, 1.0, c91, c92, c93,
		x91, y92, z92, 1.0, c91, c92, c93,

		//	right

		x92, y91, z92, 1.0, c91, c92, c93,
		x92, y91, z91, 1.0, c91, c92, c93,
		x92, y92, z92, 1.0, c91, c92, c93,
		x92, y92, z91, 1.0, c91, c92, c93,
		x92, y92, z92, 1.0, c91, c92, c93,
		x92, y91, z91, 1.0, c91, c92, c93,
		//left


		x91, y91, z92, 1.0, c91, c92, c93,
		x91, y91, z91, 1.0, c91, c92, c93,
		x91, y92, z92, 1.0, c91, c92, c93,
		x91, y92, z91, 1.0, c91, c92, c93,
		x91, y92, z92, 1.0, c91, c92, c93,
		x91, y91, z91, 1.0, c91, c92, c93,
		//	botom
		x91, y91, z91, 1.0, c91, c92, c93,
		x91, y91, z92, 1.0, c91, c92, c93,
		x92, y91, z91, 1.0, c91, c92, c93,
		x92, y91, z92, 1.0, c91, c92, c93,
		x92, y91, z91, 1.0, c91, c92, c93,
		x91, y91, z92, 1.0, c91, c92, c93,

		//	top

		x91, y92, z91, 1.0, c91, c92, c93,
		x91, y92, z92, 1.0, c91, c92, c93,
		x92, y92, z91, 1.0, c91, c92, c93,
		x92, y92, z92, 1.0, c91, c92, c93,
		x92, y92, z91, 1.0, c91, c92, c93,
		x91, y92, z92, 1.0, c91, c92, c93,



		//head end
		//body start

		// right ear
		xa1, ya1, za1, 1.0, ca1, ca2, ca3,
		xa2, ya1, za1, 1.0, ca1, ca2, ca3,
		xa1, ya2, za1, 1.0, ca1, ca2, ca3,
		xa2, ya2, za1, 1.0, ca1, ca2, ca3,
		xa2, ya1, za1, 1.0, ca1, ca2, ca3,
		xa1, ya2, za1, 1.0, ca1, ca2, ca3,

		//Face : front
		xa1, ya1, za2, 1.0, ca1, ca2, ca3,
		xa2, ya1, za2, 1.0, ca1, ca2, ca3,
		xa1, ya2, za2, 1.0, ca1, ca2, ca3,
		xa2, ya2, za2, 1.0, ca1, ca2, ca3,
		xa2, ya1, za2, 1.0, ca1, ca2, ca3,
		xa1, ya2, za2, 1.0, ca1, ca2, ca3,

		//	right

		xa2, ya1, za2, 1.0, ca1, ca2, ca3,
		xa2, ya1, za1, 1.0, ca1, ca2, ca3,
		xa2, ya2, za2, 1.0, ca1, ca2, ca3,
		xa2, ya2, za1, 1.0, ca1, ca2, ca3,
		xa2, ya2, za2, 1.0, ca1, ca2, ca3,
		xa2, ya1, za1, 1.0, ca1, ca2, ca3,

		//left


		xa1, ya1, za2, 1.0, ca1, ca2, ca3,
		xa1, ya1, za1, 1.0, ca1, ca2, ca3,
		xa1, ya2, za2, 1.0, ca1, ca2, ca3,
		xa1, ya2, za1, 1.0, ca1, ca2, ca3,
		xa1, ya2, za2, 1.0, ca1, ca2, ca3,
		xa1, ya1, za1, 1.0, ca1, ca2, ca3,

		//	bottom
		xa1, ya1, za1, 1.0, ca1, ca2, ca3,
		xa1, ya1, za2, 1.0, ca1, ca2, ca3,
		xa2, ya1, za1, 1.0, ca1, ca2, ca3,
		xa2, ya1, za2, 1.0, ca1, ca2, ca3,
		xa2, ya1, za1, 1.0, ca1, ca2, ca3,
		xa1, ya1, za2, 1.0, ca1, ca2, ca3,

		//	top

		xa1, ya2, za1, 1.0, 0, 0, 1,
		xa1, ya2, za2, 1.0, 1, 0, 0,
		xa2, ya2, za1, 1.0, 0, 1, 0,
		xa2, ya2, za2, 1.0, ca1, ca2, ca3,
		xa2, ya2, za1, 1.0, ca1, ca2, ca3,
		xa1, ya2, za2, 1.0, ca1, ca2, ca3,

		// right ear
		xb1, yb1, zb1, 1.0, cb1, cb2, cb3,
		xb2, yb1, zb1, 1.0, cb1, cb2, cb3,
		xb1, yb2, zb1, 1.0, cb1, cb2, cb3,
		xb2, yb2, zb1, 1.0, cb1, cb2, cb3,
		xb2, yb1, zb1, 1.0, cb1, cb2, cb3,
		xb1, yb2, zb1, 1.0, cb1, cb2, cb3,

		//Face : front
		xb1, yb1, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb1, zb2, 1.0, cb1, cb2, cb3,
		xb1, yb2, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb2, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb1, zb2, 1.0, cb1, cb2, cb3,
		xb1, yb2, zb2, 1.0, cb1, cb2, cb3,

		//	right

		xb2, yb1, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb1, zb1, 1.0, cb1, cb2, cb3,
		xb2, yb2, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb2, zb1, 1.0, cb1, cb2, cb3,
		xb2, yb2, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb1, zb1, 1.0, cb1, cb2, cb3,

		//left


		xb1, yb1, zb2, 1.0, cb1, cb2, cb3,
		xb1, yb1, zb1, 1.0, cb1, cb2, cb3,
		xb1, yb2, zb2, 1.0, cb1, cb2, cb3,
		xb1, yb2, zb1, 1.0, cb1, cb2, cb3,
		xb1, yb2, zb2, 1.0, cb1, cb2, cb3,
		xb1, yb1, zb1, 1.0, cb1, cb2, cb3,

		//	bottom
		xb1, yb1, zb1, 1.0, cb1, cb2, cb3,
		xb1, yb1, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb1, zb1, 1.0, cb1, cb2, cb3,
		xb2, yb1, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb1, zb1, 1.0, cb1, cb2, cb3,
		xb1, yb1, zb2, 1.0, cb1, cb2, cb3,

		//	top

		xb1, yb2, zb1, 1.0, cb1, cb2, cb3,
		xb1, yb2, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb2, zb1, 1.0, cb1, cb2, cb3,
		xb2, yb2, zb2, 1.0, cb1, cb2, cb3,
		xb2, yb2, zb1, 1.0, cb1, cb2, cb3,
		xb1, yb2, zb2, 1.0, cb1, cb2, cb3,

		//body end

		// legs start

		xc1, yc1, zc1, 1.0, cc1, cc2, cc3,
		xc2, yc1, zc1, 1.0, cc1, cc2, cc3,
		xc1, yc2, zc1, 1.0, cc1, cc2, cc3,
		xc2, yc2, zc1, 1.0, cc1, cc2, cc3,
		xc2, yc1, zc1, 1.0, cc1, cc2, cc3,
		xc1, yc2, zc1, 1.0, cc1, cc2, cc3,

		//Face : front
		xc1, yc1, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc1, zc2, 1.0, cc1, cc2, cc3,
		xc1, yc2, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc2, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc1, zc2, 1.0, cc1, cc2, cc3,
		xc1, yc2, zc2, 1.0, cc1, cc2, cc3,

		//right

		xc2, yc1, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc1, zc1, 1.0, cc1, cc2, cc3,
		xc2, yc2, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc2, zc1, 1.0, cc1, cc2, cc3,
		xc2, yc2, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc1, zc1, 1.0, cc1, cc2, cc3,

		//left


		xc1, yc1, zc2, 1.0, cc1, cc2, cc3,
		xc1, yc1, zc1, 1.0, cc1, cc2, cc3,
		xc1, yc2, zc2, 1.0, cc1, cc2, cc3,
		xc1, yc2, zc1, 1.0, cc1, cc2, cc3,
		xc1, yc2, zc2, 1.0, cc1, cc2, cc3,
		xc1, yc1, zc1, 1.0, cc1, cc2, cc3,

		//	bottom
		xc1, yc1, zc1, 1.0, cc1, cc2, cc3,
		xc1, yc1, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc1, zc1, 1.0, cc1, cc2, cc3,
		xc2, yc1, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc1, zc1, 1.0, cc1, cc2, cc3,
		xc1, yc1, zc2, 1.0, cc1, cc2, cc3,

		//	top

		xc1, yc2, zc1, 1.0, cc1, cc2, cc3,
		xc1, yc2, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc2, zc1, 1.0, cc1, cc2, cc3,
		xc2, yc2, zc2, 1.0, cc1, cc2, cc3,
		xc2, yc2, zc1, 1.0, cc1, cc2, cc3,
		xc1, yc2, zc2, 1.0, cc1, cc2, cc3,

		//Hand part
		//back
		xd1, yd1, zd1, 1.0, cd1, cd2, cd3,
		xd2, yd1, zd1, 1.0, cd1, cd2, cd3,
		xd1, yd2, zd1, 1.0, cd1, cd2, cd3,
		xd2, yd2, zd1, 1.0, cd1, cd2, cd3,
		xd2, yd1, zd1, 1.0, cd1, cd2, cd3,
		xd1, yd2, zd1, 1.0, cd1, cd2, cd3,

		//Face : front
		xd1, yd1, zd2, 1.0, cd1, cd2, cd3,
		xd2, yd1, zd2, 1.0, cd1, cd2, cd3,
		xd1, yd2, zd2, 1.0, cd1, cd2, cd3,
		xd2, yd2, zd2, 1.0, cd1, cd2, cd3,
		xd2, yd1, zd2, 1.0, cd1, cd2, cd3,
		xd1, yd2, zd2, 1.0, cd1, cd2, cd3,

		//	right

		xd2, yd1, zd2, 1.0, cd1, cd2, cd3,
		xd2, yd1, zd1, 1.0, cd1, cd2, cd3,
		xd2, yd2, zd2, 1.0, cd1, cd2, cd3,
		xd2, yd2, zd1, 1.0, cd1, cd2, cd3,
		xd2, yd2, zd2, 1.0, cd1, cd2, cd3,
		xd2, yd1, zd1, 1.0, cd1, cd2, cd3,

		//left


		xd1, yd1, zd2, 1.0, cd1, cd2, cd3,
		xd1, yd1, zd1, 1.0, cd1, cd2, cd3,
		xd1, yd2, zd2, 1.0, cd1, cd2, cd3,
		xd1, yd2, zd1, 1.0, cd1, cd2, cd3,
		xd1, yd2, zd2, 1.0, cd1, cd2, cd3,
		xd1, yd1, zd1, 1.0, cd1, cd2, cd3,

		//	bottom
		xd1, yd1, zd1, 1.0, cd1, cd2, cd3,
		xd1, yd1, zd2, 1.0, cd1, cd2, cd3,
		xd2, yd1, zd1, 1.0, cd1, cd2, cd3,
		xd2, yd1, zd2, 1.0, cd1, cd2, cd3,
		xd2, yd1, zd1, 1.0, cd1, cd2, cd3,
		xd1, yd1, zd2, 1.0, cd1, cd2, cd3,

		//	top

		xd1, yd2, zd1, 1.0, 1, 0, 0,
		xd1, yd2, zd2, 1.0, 0, 0, 1,
		xd2, yd2, zd1, 1.0, 0, 1, 0,
		xd2, yd2, zd2, 1.0, cd1, cd2, cd3,
		xd2, yd2, zd1, 1.0, cd1, cd2, cd3,
		xd1, yd2, zd2, 1.0, cd1, cd2, cd3,

	]);




}

function calc_lookat(){

	if (w_is_pushed){
		x_location = x_location + 0.1
	}

	if (s_is_pushed){

		x_location = x_location - 0.1
	}


	return [x_location, y_location, z_location, x_refrence, y_refrence, z_refrence, x_up, y_up, z_up]
}

function draw_two(gl, n, currentAngle, modelMatrix, u_ModelMatrix){

	

	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

	var aspect = gl.canvas.width / 2 /gl.canvas.height
	var angle = 35.0
	var lookat_val = calc_lookat();
	view.setLookAt(
		lookat_val[0],
		lookat_val[1],
		lookat_val[2],
		lookat_val[3],
		lookat_val[4],
		lookat_val[5],
		lookat_val[6],
		lookat_val[7],
		lookat_val[8]);
	
	
	projection.setPerspective(angle, 1 ,aspect, 100.0);
	
	
	gl.uniformMatrix4fv(costum_Projmatrix, false, projection.elements);
	gl.uniformMatrix4fv(costum_Viewmatrix, false, view.elements);

	gl.viewport(0, 0, gl.canvas.width / 2, gl.canvas.height);
	drawAll(gl, n, currentAngle, modelMatrix, u_ModelMatrix);    // Draw shapes
	
	var width = 3.0
	var height = width / aspect
	projection.setOrtho(-width/2,width/2, -height/2, height/2, 0.0, 100.0)

	gl.uniformMatrix4fv(costum_Projmatrix, false, projection.elements);
	gl.viewport(gl.canvas.width / 2, 0, gl.canvas.width / 2, gl.canvas.height);
	drawAll(gl, n, currentAngle, modelMatrix, u_ModelMatrix);   // Draw shapes


}

function drawAll(gl, n, currentAngle, modelMatrix, u_ModelMatrix) {
	//==============================================================================
	// Clear <canvas>  colors AND the depth buffer
	
	modelMatrix.setIdentity();    // DEFINE 'world-space' coords.

	//draw world space coordinates
	pushMatrix(modelMatrix);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	
	gl.drawArrays(gl.LINES, 								
		axesStart / floatsPerVertex,	
		axes_vars.length / floatsPerVertex);


	popMatrix(modelMatrix);



	pushMatrix(modelMatrix);  // SAVE world drawing coords.
	//---------Draw Ground Plane, without spinning.
	// position it.
	modelMatrix.translate(0.4, -0.4, 0.0);
	modelMatrix.scale(0.1, 0.1, 0.1);				// shrink by 10X:

	// Drawing:
	// Pass our current matrix to the vertex shaders:
	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw just the ground-plane's vertices
	gl.drawArrays(gl.LINES, 								// use this drawing primitive, and
		gndStart / floatsPerVertex,	// start at this vertex number, and
		gndVerts.length / floatsPerVertex);	// draw this many vertices.
	modelMatrix = popMatrix();  // RESTORE 'world' drawing coords.

	//make figures stand up

	
	modelMatrix.rotate(90,1,0,0);
	modelMatrix.translate(0,1,0)
 	//new shape
	 pushMatrix(modelMatrix);     // SAVE world coord system;
	 //-------Draw Spinning Cylinder:
	 
	 modelMatrix.translate(-3, -0.4, 0.0);  // 'set' means DISCARD old matrix,
	 // (drawing axes centered in CVV), and then make new
	 // drawing axes moved to the lower-left corner of CVV. 
	 modelMatrix.scale(0.2, 0.2, 0.2);
	 // if you DON'T scale, cyl goes outside the CVV; clipped!
	 modelMatrix.rotate(270, 1, 0, 0);  // spin around y axis.
	 // Drawing:
	 // Pass our current matrix to the vertex shaders:
	 gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	 // Draw the cylinder's vertices, and no other vertices:
	 gl.drawArrays(gl.TRIANGLE_STRIP,				// use this drawing primitive, and
		 cylStart / floatsPerVertex, // start at this vertex number, and
		 cylVerts.length / floatsPerVertex);	// draw this many vertices.
	 modelMatrix = popMatrix();  // RESTORE 'world' drawing coords.
	 //===========================================================
	 //  
	 pushMatrix(modelMatrix);  // SAVE world drawing coords.
	 //--------Draw Spinning torus
	 modelMatrix.translate(-3, -0.1, 0);	// 'set' means DISCARD old matrix,
 
	 modelMatrix.scale(0.3, 0.2, 0.3);
	 // Make it smaller:
	 modelMatrix.rotate(45 + g_angle02, 1, 0, 0);  // Spin on YZ axis
	 // Drawing:		
	 // Pass our current matrix to the vertex shaders:
	 gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	 // Draw just the torus's vertices
	 gl.drawArrays(gl.TRIANGLE_STRIP, 				// use this drawing primitive, and
		 torStart / floatsPerVertex,	// start at this vertex number, and
		 torVerts.length / floatsPerVertex);	// draw this many vertices.
	 modelMatrix = popMatrix();  // RESTORE 'world' drawing coords.
	 //===========================================================
	









	//===========================================================
	//old shapes
	modelMatrix.setTranslate(-0.25, 0, 0.3);

	// convert to left-handed coord sys


	modelMatrix.rotate(90,1,0,0);



	var dist = Math.sqrt(g_xMdragTot * g_xMdragTot + g_yMdragTot * g_yMdragTot);

	modelMatrix.rotate(dist * 120.0, -g_yMdragTot + 0.0001, g_xMdragTot + 0.0001, 0.0);

	modelMatrix.scale(0.5, 0.5, 0.5);

	pushMatrix(modelMatrix);
	modelMatrix.rotate(-g_angle05, 0, 1, 0);
	modelMatrix.translate(0.2, 0, 0);
	pushMatrix(modelMatrix);



	//Body
	modelMatrix.rotate(-g_angle02, 1, 0, 0);
	modelMatrix.translate(0, 0, -0.1);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 180, 36);


	//head

	pushMatrix(modelMatrix);

	modelMatrix.scale(0.75, 0.75, 0.75);

	modelMatrix.translate(0.125, 1.075, 0);




	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 0, 108);
	modelMatrix = popMatrix();

	//right arm

	//	g_modelMatrix = popMatrix(g_modelMatrix);
	pushMatrix(modelMatrix);

	modelMatrix.scale(0.5, 0.5, 0.5);
	modelMatrix.translate(0, 1.55, 0);



	modelMatrix.rotate(g_angle01 + 180, 0, 0, 1);
	modelMatrix.translate(-0.1, -0.5, 0);



	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 108, 72);

	//draw stick
	modelMatrix.scale(0.75, 0.75, 0.75);

	modelMatrix.translate(0, -0.3, 0.2);
	modelMatrix.rotate( g_angle02 + 180, 1, 0, 0);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 216, 36);


	modelMatrix = popMatrix();
	//left arm

	pushMatrix(modelMatrix);

	modelMatrix.scale(0.5, 0.5, 0.5);
	modelMatrix.translate(0.75, 1.50, 0.25);

	modelMatrix.rotate(g_angle01, 0, 0, 1);
	modelMatrix.translate(-0.1, -0.5, 0);


	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 108, 72);
	

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	
	gl.drawArrays(gl.LINES, 								
		axesStart / floatsPerVertex,	
		axes_vars.length / floatsPerVertex);

	modelMatrix = popMatrix();




	modelMatrix = popMatrix();


	//right leg
	pushMatrix(modelMatrix);

	modelMatrix.scale(0.75, 0.75, 0.75);

	modelMatrix.translate(0.01, -0.7, 0);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 216, 72);

	modelMatrix = popMatrix();
	//left leg
	pushMatrix(modelMatrix);

	modelMatrix.scale(0.75, 0.75, 0.75);

	modelMatrix.translate(0.32, -0.7, 0);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 216, 72);


	modelMatrix = popMatrix();
	modelMatrix = popMatrix();

	//cowboy end

	//pig start

	modelMatrix.rotate(g_angle04, 0, 1, 0);
	modelMatrix.translate(0.6, -0.2, 0);

	pushMatrix(modelMatrix);

	//body


	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 396, 36);

	modelMatrix = popMatrix();

	//body end
	//face start
	pushMatrix(modelMatrix);
	modelMatrix.translate(0.05, 0.2, 0);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 288, 108);
	modelMatrix = popMatrix();

	//forder right leg

	pushMatrix(modelMatrix);

	modelMatrix.scale(0.75, 0.75, 0.75);

	modelMatrix.translate(0, 0, 0.1);

	modelMatrix.rotate(g_angle03, 1, 0, 0);
	modelMatrix.translate(0, -0.35, -0.1);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2wsw triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 432, 72);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	
	gl.drawArrays(gl.LINES, 								
		axesStart / floatsPerVertex,	
		axes_vars.length / floatsPerVertex);




	modelMatrix = popMatrix();

	//forder left leg

	pushMatrix(modelMatrix);

	modelMatrix.scale(0.75, 0.75, 0.75);

	modelMatrix.translate(0.33, 0, 0.1);

	modelMatrix.rotate(g_angle03, 1, 0, 0);
	modelMatrix.translate(0, -0.35, -0.1);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 432, 72);
	modelMatrix = popMatrix();

	//back right leg

	pushMatrix(modelMatrix);

	modelMatrix.scale(0.75, 0.75, 0.75);

	modelMatrix.translate(0, 0, 0.7);

	modelMatrix.rotate(-g_angle03, 1, 0, 0);
	modelMatrix.translate(0, -0.35, -0.1);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 432, 72);
	modelMatrix = popMatrix();

	//back left leg
	pushMatrix(modelMatrix);

	modelMatrix.scale(0.75, 0.75, 0.75);

	modelMatrix.translate(0.33, 0, 0.7);

	modelMatrix.rotate(-g_angle03, 1, 0, 0);
	modelMatrix.translate(0, -0.35, -0.1);

	gl.uniformMatrix4fv(u_ModelMatrix, false, modelMatrix.elements);
	// Draw only the last 2 triangles: start at vertex 6, draw 6 vertices
	gl.drawArrays(gl.TRIANGLES, 432, 72);
	modelMatrix = popMatrix();	

}

// Last time that this function was called:  (used for animation timing)
var g_last = Date.now();

function animate(angle) {
	//==============================================================================
	// Calculate the elapsed time
	var now = Date.now();
	var elapsed = now - g_last;
	g_last = now;
	// Update the current rotation angle (adjusted by the elapsed time)
	//  limit the angle to move smoothly between +20 and -85 degrees:
	//  if(angle >  120.0 && ANGLE_STEP > 0) ANGLE_STEP = -ANGLE_STEP;
	//  if(angle < -120.0 && ANGLE_STEP < 0) ANGLE_STEP = -ANGLE_STEP;


	//control of cowboy

	



	if (g_isRun == true) {
		g_angle05 = g_angle05 + (g_angle05Rate * elapsed) / 1000.0;
	}

	g_angle05 %= 360;
	if (g_isRun == true) {
		g_angle01 = g_angle01 + (g_angle01Rate * elapsed) / 1000.0;
	}

	if (g_angle01 > 160.0 && g_angle01Rate > 0) g_angle01Rate = -g_angle01Rate
	if (g_angle01 < 10.0 && g_angle01Rate < 0) g_angle01Rate = -g_angle01Rate

	if (ArrowUP_is_pushed){

		g_angle01Rate = g_angle01Rate * 1.02;
		g_angle02Rate = g_angle02Rate * 1.02;
	}
	if (ArrowDOWN_is_pushed){

		g_angle01Rate = g_angle01Rate * 0.98;
		g_angle02Rate = g_angle02Rate * 0.98;
	}


	if (g_isRun == true) {

		g_angle02 = g_angle02 + (g_angle02Rate * elapsed) / 1000.0;
	}

	if (g_angle02 > 180.0) g_angle02 = g_angle02 - 360.0;
	if (g_angle02 < -180.0) g_angle02 = g_angle02 + 360.0;

	if (g_angle02 > 60.0 && g_angle02Rate > 0) g_angle02Rate *= -1.0;
	if (g_angle02 < 0.0 && g_angle02Rate < 0) g_angle02Rate *= -1.0;

	var eps = 5



	if ((w_is_pushed || s_is_pushed) || !((g_angle03 < eps) && ((g_angle03 > -eps)))) {
		g_angle03 = g_angle03 + (g_angle03Rate * elapsed) / 1000.0;
		if (g_angle03 > 100.0 && g_angle03Rate > 0) g_angle03Rate = -g_angle03Rate
		if (g_angle03 < -10.0 && g_angle03Rate < 0) g_angle03Rate = -g_angle03Rate
	}

	if ((w_is_pushed || s_is_pushed) || !((g_angle03 < eps) && ((g_angle03 > -eps)))) {

		if (w_is_pushed || w_last_pushed) {
			g_angle04 = g_angle04 + (g_angle04Rate * elapsed) / 1000.0;
			w_last_pushed = true;
		}
		else {
			g_angle04 = g_angle04 - (g_angle04Rate * elapsed) / 1000.0;
			w_last_pushed = false
		}

	}
	var newAngle = angle + (ANGLE_STEP * elapsed) / 1000.0;
	return newAngle %= 360;
}

//==================HTML Button Callbacks
function nextShape() {
	shapeNum += 1;
	if (shapeNum >= shapeMax) shapeNum = 0;
}

function spinDown() {
	ANGLE_STEP -= 25;
}

function spinUp() {
	ANGLE_STEP += 25;
}

function runStop() {
	if (ANGLE_STEP * ANGLE_STEP > 1) {
		myTmp = ANGLE_STEP;
		ANGLE_STEP = 0;
	}
	else {
		ANGLE_STEP = myTmp;
	}
}

//===================Mouse and Keyboard event-handling Callbacks

function myMouseDown(ev) {
	//==============================================================================
	// Called when user PRESSES down any mouse button;
	// 									(Which button?    console.log('ev.button='+ev.button);   )
	// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
	//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)  

	// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
	var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
	var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
	var yp = g_canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
	//  console.log('myMouseDown(pixel coords): xp,yp=\t',xp,',\t',yp);

	// Convert to Canonical View Volume (CVV) coordinates too:
	var x = (xp - g_canvas.width / 2) / 		// move origin to center of canvas and
		(g_canvas.width / 2);			// normalize canvas to -1 <= x < +1,
	var y = (yp - g_canvas.height / 2) /		//										 -1 <= y < +1.
		(g_canvas.height / 2);
	//	console.log('myMouseDown(CVV coords  ):  x, y=\t',x,',\t',y);

	g_isDrag = true;											// set our mouse-dragging flag
	g_xMclik = x;													// record where mouse-dragging began
	g_yMclik = y;
	// report on webpage
	document.getElementById('MouseAtResult').innerHTML =
		'Mouse At: ' + x.toFixed(g_digits) + ', ' + y.toFixed(g_digits);
};


function myMouseMove(ev) {
	//==============================================================================
	// Called when user MOVES the mouse with a button already pressed down.
	// 									(Which button?   console.log('ev.button='+ev.button);    )
	// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
	//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)  

	if (g_isDrag == false) return;				// IGNORE all mouse-moves except 'dragging'

	// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
	var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
	var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
	var yp = g_canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
	//  console.log('myMouseMove(pixel coords): xp,yp=\t',xp,',\t',yp);

	// Convert to Canonical View Volume (CVV) coordinates too:
	var x = (xp - g_canvas.width / 2) / 		// move origin to center of canvas and
		(g_canvas.width / 2);		// normalize canvas to -1 <= x < +1,
	var y = (yp - g_canvas.height / 2) /		//									-1 <= y < +1.
		(g_canvas.height / 2);
	//	console.log('myMouseMove(CVV coords  ):  x, y=\t',x,',\t',y);

	// find how far we dragged the mouse:
	g_xMdragTot += (x - g_xMclik);			// Accumulate change-in-mouse-position,&
	g_yMdragTot += (y - g_yMclik);
	// Report new mouse position & how far we moved on webpage:
	document.getElementById('MouseAtResult').innerHTML =
		'Mouse At: ' + x.toFixed(g_digits) + ', ' + y.toFixed(g_digits);
	document.getElementById('MouseDragResult').innerHTML =
		'Mouse Drag: ' + (x - g_xMclik).toFixed(g_digits) + ', '
		+ (y - g_yMclik).toFixed(g_digits);

	g_xMclik = x;											// Make next drag-measurement from here.
	g_yMclik = y;
};

function myMouseUp(ev) {
	//==============================================================================
	// Called when user RELEASES mouse button pressed previously.
	// 									(Which button?   console.log('ev.button='+ev.button);    )
	// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
	//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!)  

	// Create right-handed 'pixel' coords with origin at WebGL canvas LOWER left;
	var rect = ev.target.getBoundingClientRect();	// get canvas corners in pixels
	var xp = ev.clientX - rect.left;									// x==0 at canvas left edge
	var yp = g_canvas.height - (ev.clientY - rect.top);	// y==0 at canvas bottom edge
	//  console.log('myMouseUp  (pixel coords):\n\t xp,yp=\t',xp,',\t',yp);

	// Convert to Canonical View Volume (CVV) coordinates too:
	var x = (xp - g_canvas.width / 2) / 		// move origin to center of canvas and
		(g_canvas.width / 2);			// normalize canvas to -1 <= x < +1,
	var y = (yp - g_canvas.height / 2) /		//										 -1 <= y < +1.
		(g_canvas.height / 2);
	console.log('myMouseUp  (CVV coords  ):\n\t x, y=\t', x, ',\t', y);

	g_isDrag = false;											// CLEAR our mouse-dragging flag, and
	// accumulate any final bit of mouse-dragging we did:
	g_xMdragTot += (x - g_xMclik);
	g_yMdragTot += (y - g_yMclik);
	// Report new mouse position:
	document.getElementById('MouseAtResult').innerHTML =
		'Mouse At: ' + x.toFixed(g_digits) + ', ' + y.toFixed(g_digits);
	console.log('myMouseUp: g_xMdragTot,g_yMdragTot =',
		g_xMdragTot.toFixed(g_digits), ',\t', g_yMdragTot.toFixed(g_digits));
};

function myMouseClick(ev) {
	//=============================================================================
	// Called when user completes a mouse-button single-click event 
	// (e.g. mouse-button pressed down, then released)
	// 									   
	//    WHICH button? try:  console.log('ev.button='+ev.button); 
	// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
	//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!) 
	//    See myMouseUp(), myMouseDown() for conversions to  CVV coordinates.

	// STUB
	console.log("myMouseClick() on button: ", ev.button);
}

function myMouseDblClick(ev) {
	//=============================================================================
	// Called when user completes a mouse-button double-click event 
	// 									   
	//    WHICH button? try:  console.log('ev.button='+ev.button); 
	// 		ev.clientX, ev.clientY == mouse pointer location, but measured in webpage 
	//		pixels: left-handed coords; UPPER left origin; Y increases DOWNWARDS (!) 
	//    See myMouseUp(), myMouseDown() for conversions to  CVV coordinates.

	// STUB
	console.log("myMouse-DOUBLE-Click() on button: ", ev.button);
}

function myKeyDown(kev) {
	//===============================================================================
	// Called when user presses down ANY key on the keyboard;
	//
	// For a light, easy explanation of keyboard events in JavaScript,
	// see:    http://www.kirupa.com/html5/keyboard_events_in_javascript.htm
	// For a thorough explanation of a mess of JavaScript keyboard event handling,
	// see:    http://javascript.info/tutorial/keyboard-events
	//
	// NOTE: Mozilla deprecated the 'keypress' event entirely, and in the
	//        'keydown' event deprecated several read-only properties I used
	//        previously, including kev.charCode, kev.keyCode. 
	//        Revised 2/2019:  use kev.key and kev.code instead.
	//
	// Report EVERYTHING in console:
	console.log("--kev.code:", kev.code, "\t\t--kev.key:", kev.key,
		"\n--kev.ctrlKey:", kev.ctrlKey, "\t--kev.shiftKey:", kev.shiftKey,
		"\n--kev.altKey:", kev.altKey, "\t--kev.metaKey:", kev.metaKey);

	// and report EVERYTHING on webpage:

	// key details:


	switch (kev.code) {
		case "KeyP":
			console.log("Pause/unPause!\n");                // print on console,
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown() found p/P key. Pause/unPause!';   // print on webpage
			if (g_isRun == true) {
				g_isRun = false;    // STOP animation
			}
			else {
				g_isRun = true;     // RESTART animation
				tick();
			}
			break;
		//------------------WASD navigation-----------------
		case "KeyA":
			console.log("a/A key: Strafe LEFT!\n");
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown() found a/A key. Strafe LEFT!';
			break;
		case "KeyD":
			console.log("d/D key: Strafe RIGHT!\n");
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown() found d/D key. Strafe RIGHT!';
			break;
		case "KeyS":
			console.log("s/S key: Move BACK!\n");
			s_is_pushed = true;
			w_last_pushed = false;
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown() found s/Sa key. Move BACK.';
			break;
		case "KeyW":
			w_is_pushed = true;
			w_last_pushed = true;
			console.log("w/W key: Move FWD!\n");
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown() found w/W key. Move FWD!';
			break;
		//----------------Arrow keys------------------------
		case "ArrowLeft":
			console.log(' left-arrow.');
			// and print on webpage in the <div> element with id='Result':
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown(): Left Arrow=' + kev.keyCode;
			break;
		case "ArrowRight":
			console.log('right-arrow.');
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown():Right Arrow:keyCode=' + kev.keyCode;
			break;
		case "ArrowUp":
			ArrowUP_is_pushed = true;
			console.log('   up-arrow.');
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown():   Up Arrow:keyCode=' + kev.keyCode;
			break;
		case "ArrowDown":
			ArrowDOWN_is_pushed = true;
			console.log(' down-arrow.');
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown(): Down Arrow:keyCode=' + kev.keyCode;
			break;
		default:
			console.log("UNUSED!");
			document.getElementById('KeyDownResult').innerHTML =
				'myKeyDown(): UNUSED!';
			break;
	}
}

function myKeyUp(kev) {
	//===============================================================================
	// Called when user releases ANY key on the keyboard; captures scancodes well
	s_is_pushed = false;
	w_is_pushed = false;
	ArrowUP_is_pushed = false;
	ArrowDOWN_is_pushed = false;
	console.log('myKeyUp()--keyCode=' + kev.keyCode + ' released.');
}

