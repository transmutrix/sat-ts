// Charlie Dean - https://github.com/transmutrix/sat-ts
// This is a TypeScript rewrite of sat-js (https://github.com/jriecken/sat-js) by Jim Riecken.
//
// Released under the MIT License - https://github.com/jriecken/sat-js
//
// A simple library for determining intersections of circles and
// polygons using the Separating Axis Theorem.

/**
 * Represents a vector in two dimensions with `x` and `y` properties.
 */
 export class Vector {
    /**
     * Create a new Vector, optionally passing in the `x` and `y` coordinates. If
     * a coordinate is not specified, it will be set to `0`
     * @param {?number=} x The x position.
     * @param {?number=} y The y position.
     * @constructor
     */
    constructor(
        public x: number = 0,
        public y: number = 0
    ) {}

    /**
     * Copy the values of another Vector into this one.
     * @param other The other Vector.
     * @returns this for chaining.
     */
    copy(other: Vector) {
        this.x = other.x;
        this.y = other.y;
        return this;
    }

    /**
     * Create a new vector with the same coordinates as this on.
     * @returns The new cloned vector
     */
    clone() {
        return new Vector(this.x, this.y);
    }
    
    /**
     * Change this vector to be perpendicular to what it was before. (Effectively
     * rotates it 90 degrees in a clockwise direction)
     * @returns This for chaining.
     */
    perp() {
        const {x} = this;
        this.x = this.y;
        this.y = -x;
        return this;
    }

    /**
     * Rotate this vector (counter-clockwise) by the specified angle (in radians).
     * @param angle The angle to rotate (in radians)
     * @returns This for chaining.
     */
    rotate(radians: number) {
        const {x, y} = this;
        this.x = x * Math.cos(radians) - y * Math.sin(radians);
        this.y = x * Math.sin(radians) + y * Math.cos(radians);
        return this;
    }

    /**
     * Reverse this vector.
     * @returns This for chaining.
     */
    reverse() {
        this.x = -this.x;
        this.y = -this.y;
        return this;
    }

    /**
     * Normalize this vector.  (make it have length of `1`)
     * @returns This for chaining.
     */
    normalize() {
        const d = this.len();
        if (d > 0) {
            this.x /= d;
            this.y /= d;
        }
        return this;
    }

    /**
     * Add another vector to this one.
     * @param other The other Vector.
     * @returns This for chaining.
     */
    add(other: Vector) {
        this.x += other.x;
        this.y += other.y;
        return this;
    }

    /**
     * Subtract another vector from this one.
     * @param {Vector} other The other Vector.
     * @returns {Vector} This for chaiing.
     */
    sub(other: Vector) {
        this.x -= other.x;
        this.y -= other.y;
        return this;
    }

    /**
     * Scale this vector. An independent scaling factor can be provided
     * for each axis, or a single scaling factor that will scale both `x` and `y`.
     * @param x The scaling factor in the x direction.
     * @param y The scaling factor in the y direction. If this is not specified, the x scaling factor will be used.
     * @returns This for chaining.
     */
    scale(x: number, y = x) {
        this.x *= x;
        this.y *= y;
        return this;
    }

    /**
     * Project this vector on to another vector.
     * @param other The vector to project onto.
     * @returns This for chaining.
     */
    project(other: Vector) {
        const amt = this.dot(other) / other.len2();
        this.x = amt * other.x;
        this.y = amt * other.y;
        return this;
    }

    /**
     * Project this vector onto a vector of unit length. This is slightly more efficient
     * than `project` when dealing with unit vectors.
     * @param other The unit vector to project onto.
     * @returns This for chaining.
     */
    projectN(other: Vector) {
        const amt = this.dot(other);
        this.x = amt * other.x;
        this.y = amt * other.y;
        return this;
    }

    /**
     * Reflect this vector on an arbitrary axis.
     * @param axis The vector representing the axis.
     * @returns This for chaining.
     */
    reflect(axis: Vector) {
        const {x, y} = this;
        this.project(axis).scale(2);
        this.x -= x;
        this.y -= y;
        return this;
    }

    /**
     * Reflect this vector on an arbitrary axis (represented by a unit vector). This is
     * slightly more efficient than `reflect` when dealing with an axis that is a unit vector.
     * @param axis The unit vector representing the axis.
     * @returns This for chaining.
     */
    reflectN(axis: Vector) {
        const {x, y} = this;
        this.projectN(axis).scale(2);
        this.x -= x;
        this.y -= y;
        return this;
    }

    /**
     * Get the dot product of this vector and another.
     * @param other The vector to dot this one against.
     * @returns The dot product.
     */
    dot(other: Vector) {
        return this.x * other.x + this.y * other.y;
    }

    /**
     * Get the squared length of this vector.
     * @returns The length^2 of this vector.
     */
    len2() {
        return this.dot(this);
    }

    /**
     * Get the length of this vector.
     * @returns The length of this vector.
     */
    len() {
        return Math.sqrt(this.len2());
    }
}

/**
 * Represents an axis-aligned box, with a width and height.
 */
export class Box {
    /**
     * Create a new box, with the specified position, width, and height. If no position
     * is given, the position will be `(0,0)`. If no width or height are given, they will
     * be set to `0`.
     * @param pos A vector representing the bottom-left of the box (i.e. the smallest x and smallest y value).
     * @param w The width of the box.
     * @param h The height of the box.
     * @constructor
     */
    constructor(
        public pos = new Vector(),
        public w = 0,
        public h = 0
    ) {}

    /**
     * Returns a polygon whose edges are the same as this box.
     * @returns A new Polygon that represents this box.
     */
    toPolygon() {
        const pos = this.pos;
        const w = this.w;
        const h = this.h;
        return new Polygon(new Vector(pos.x, pos.y), [
            new Vector(0, 0), new Vector(w, 0),
            new Vector(w, h), new Vector(0, h)
        ]);
    }
}

/**
 * Represents a circle with a position and a radius.
 */
export class Circle {
    offset = new Vector();
    
    /**
     * Create a new circle, optionally passing in a position and/or radius. If no position
     * is given, the circle will be at `(0,0)`. If no radius is provided, the circle will
     * have a radius of `0`.
     * @param pos A vector representing the position of the center of the circle
     * @param r The radius of the circle
     * @constructor
     */
    constructor(
        public pos = new Vector(),
        public r = 0
    ) {}

    /**
     * Compute the axis-aligned bounding box (AABB) of this Circle.
     * **Note: Returns a _new_ `Box` each time you call this.**
     * @returns The AABB
     */
    getAABBAsBox() {
        const r = this.r;
        const corner = this.pos.clone().add(this.offset).sub(new Vector(r, r));
        return new Box(corner, r * 2, r * 2);
    }

    /**
     * Compute the axis-aligned bounding box (AABB) of this Circle.
     * **Note: Returns a _new_ `Polygon` each time you call this.**
     * @returns The AABB
     */
    getAABB() {
        return this.getAABBAsBox().toPolygon();
    }

    /**
     * Set the current offset to apply to the radius.
     * @param offset The new offset vector.
     * @returns This for chaining.
     */
    setOffset(offset: Vector) {
        this.offset = offset;
        return this;
    }
}

/**
 * Represents a *convex* polygon with any number of points (specified in counter-clockwise order)
 * Note: Do _not_ manually change the `points`, `angle`, or `offset` properties. Use the
 * provided setters. Otherwise the calculated properties will not be updated correctly.
 * `pos` can be changed directly.
 */
export class Polygon {
    points: Vector[] = [];
    angle = 0;
    offset = new Vector();

    calcPoints: Vector[] = [];
    edges: Vector[] = [];
    normals: Vector[] = [];

    /**
     * Create a new polygon, passing in a position vector, and an array of points (represented
     * by vectors relative to the position vector). If no position is passed in, the position
     * of the polygon will be `(0,0)`.
     * @param pos A vector representing the origin of the polygon. (all other points are relative to this one)
     * @param points An array of vectors representing the points in the polygon, in counter-clockwise order.
     * @constructor
     */
    constructor(
        public pos = new Vector(),
        points = new Array<Vector>()
    ) {
        this.setPoints(points);
    }

    /**
     * Set the points of the polygon. Any consecutive duplicate points will be combined.
     *
     * Note: The points are counter-clockwise *with respect to the coordinate system*.
     * If you directly draw the points on a screen that has the origin at the top-left corner
     * it will _appear_ visually that the points are being specified clockwise. This is just
     * because of the inversion of the Y-axis when being displayed.
     * @param points An array of vectors representing the points in the polygon, in counter-clockwise order.
     * @returns This for chaining.
     */
    setPoints(points: Vector[]) {
        // Only re-allocate if this is a new polygon or the number of points has changed.
        const lengthChanged = !this.points || this.points.length !== points.length;
        if (lengthChanged) {
            let i;
            const calcPoints: Vector[] = this.calcPoints = [];
            const edges: Vector[] = this.edges = [];
            const normals: Vector[] = this.normals = [];
            // Allocate the vector arrays for the calculated properties
            for (i = 0; i < points.length; i++) {
                // Remove consecutive duplicate points
                const p1 = points[i];
                const p2 = i < points.length - 1 ? points[i + 1] : points[0];
                if (p1 !== p2 && p1.x === p2.x && p1.y === p2.y) {
                    points.splice(i, 1);
                    i -= 1;
                    continue;
                }
                calcPoints.push(new Vector());
                edges.push(new Vector());
                normals.push(new Vector());
            }
        }
        this.points = points;
        this._recalc();
        return this;
    }

    /**
     * Set the current rotation angle of the polygon.
     * @param angle The current rotation angle (in radians).
     * @returns This for chaining.
     */
    setAngle(radians: number) {
        this.angle = radians;
        this._recalc();
        return this;
    }

    /**
     * Set the current offset to apply to the `points` before applying the `angle` rotation.
     * @param offset The new offset vector.
     * @returns This for chaining.
     */
    setOffset(offset: Vector) {
        this.offset = offset;
        this._recalc();
        return this;
    }

    /**
     * Rotates this polygon counter-clockwise around the origin of *its local coordinate system* (i.e. `pos`).
     * **Note: This changes the **original** points (so any `angle` will be applied on top of this rotation).**
     * @param angle The angle to rotate (in radians)
     * @returns This for chaining.
     */
    rotate(radians: number) {
        const points = this.points;
        const len = points.length;
        for (let i = 0; i < len; i++) {
            points[i].rotate(radians);
        }
        this._recalc();
        return this;
    }

    /**
     * Translates the points of this polygon by a specified amount relative to the origin of *its own coordinate
     * system* (i.e. `pos`).
     * This is most useful to change the "center point" of a polygon. If you just want to move the whole polygon, change
     * the coordinates of `pos`.
     * Note: This changes the **original** points (so any `offset` will be applied on top of this translation)
     * @param x The horizontal amount to translate.
     * @param y The vertical amount to translate.
     * @returns This for chaining.
     */
    translate(x: number, y: number) {
        const points = this.points;
        const len = points.length;
        for (let i = 0; i < len; i++) {
            points[i].x += x;
            points[i].y += y;
        }
        this._recalc();
        return this;
    }

    /**
     * Computes the calculated collision polygon. Applies the `angle` and `offset` to the original points then recalculates the
     * edges and normals of the collision polygon.
     * @returns This for chaining.
     */
    private _recalc() {
        // Calculated points - this is what is used for underlying collisions and takes into account
        // the angle/offset set on the polygon.
        const calcPoints = this.calcPoints;
        // The edges here are the direction of the `n`th edge of the polygon, relative to
        // the `n`th point. If you want to draw a given edge from the edge value, you must
        // first translate to the position of the starting point.
        const edges = this.edges;
        // The normals here are the direction of the normal for the `n`th edge of the polygon, relative
        // to the position of the `n`th point. If you want to draw an edge normal, you must first
        // translate to the position of the starting point.
        const normals = this.normals;
        // Copy the original points array and apply the offset/angle
        const points = this.points;
        const offset = this.offset;
        const angle = this.angle;
        const len = points.length;
        let i;
        for (i = 0; i < len; i++) {
            const calcPoint = calcPoints[i].copy(points[i]);
            calcPoint.x += offset.x;
            calcPoint.y += offset.y;
            if (angle !== 0) {
                calcPoint.rotate(angle);
            }
        }
        // Calculate the edges/normals
        for (i = 0; i < len; i++) {
            const p1 = calcPoints[i];
            const p2 = i < len - 1 ? calcPoints[i + 1] : calcPoints[0];
            const e = edges[i].copy(p2).sub(p1);
            normals[i].copy(e).perp().normalize();
        }
        return this;
    }

    /**
     * Compute the axis-aligned bounding box. Any current state
     * (translations/rotations) will be applied before constructing the AABB.
     * **Note: Returns a _new_ `Box` each time you call this.**
     * @returns The AABB
     */
    getAABBAsBox() {
        const points = this.calcPoints;
        const len = points.length;
        let xMin = points[0].x;
        let yMin = points[0].y;
        let xMax = points[0].x;
        let yMax = points[0].y;
        for (let i = 1; i < len; i++) {
            const point = points[i];
            if (point.x < xMin) {
                xMin = point.x;
            }
            else if (point.x > xMax) {
                xMax = point.x;
            }
            if (point.y < yMin) {
                yMin = point.y;
            }
            else if (point.y > yMax) {
                yMax = point.y;
            }
        }
        return new Box(this.pos.clone().add(new Vector(xMin, yMin)), xMax - xMin, yMax - yMin);
    }

    /**
     * Compute the axis-aligned bounding box. Any current state
     * (translations/rotations) will be applied before constructing the AABB.
     * **Note: Returns a _new_ `Polygon` each time you call this.**
     * @returns The AABB
     */
    getAABB() {
        return this.getAABBAsBox().toPolygon();
    }

    /**
     * Compute the centroid (geometric center) of the polygon. Any current state
     * (translations/rotations) will be applied before computing the centroid.
     * See https://en.wikipedia.org/wiki/Centroid#Centroid_of_a_polygon
     * **Note: Returns a _new_ `Vector` each time you call this.**
     * @returns A Vector that contains the coordinates of the Centroid.
     */
    getCentroid() {
        const points = this.calcPoints;
        const len = points.length;
        let cx = 0;
        let cy = 0;
        let ar = 0;
        for (let i = 0; i < len; i++) {
            const p1 = points[i];
            const p2 = i === len - 1 ? points[0] : points[i + 1]; // Loop around if last point
            const a = p1.x * p2.y - p2.x * p1.y;
            cx += (p1.x + p2.x) * a;
            cy += (p1.y + p2.y) * a;
            ar += a;
        }
        ar = ar * 3; // we want 1 / 6 the area and we currently have 2*area
        cx = cx / ar;
        cy = cy / ar;
        return new Vector(cx, cy);
    }
}

/**
 * An object representing the result of an intersection. Contains:
 * - The two objects participating in the intersection
 * - The vector representing the minimum change necessary to extract the first object
 *   from the second one (as well as a unit vector in that direction and the magnitude
 *   of the overlap)
 * - Whether the first object is entirely inside the second, and vice versa.
 */
export class Response {
    a?: Circle|Box|Polygon;
    b?: Circle|Box|Polygon;
    overlapN = new Vector();
    overlapV = new Vector();
    overlap = Number.MAX_VALUE;
    aInB = true;
    bInA = true;

    /**
     * @constructor
     */
    constructor() {
        this.clear();
    }

    /**
     * Set some values of the response back to their defaults.  Call this between tests if
     * you are going to reuse a single Response object for multiple intersection tests (recommented
     * as it will avoid allcating extra memory)
     * @returns This for chaining
     */
    clear() {
        this.aInB = true;
        this.bInA = true;
        this.overlap = Number.MAX_VALUE;
        return this;
    }
}

/// Object Pools

/**
 * A pool of `Vector` objects that are used in calculations to avoid
 * allocating memory.
 */
const T_VECTORS: Vector[] = [];
for (let i = 0; i < 10; i++) { T_VECTORS.push(new Vector()); }

/**
 * A pool of arrays of numbers used in calculations to avoid allocating memory.
 */
const T_ARRAYS: number[][] = [];
for (let i = 0; i < 5; i++) { T_ARRAYS.push([]); }

/**
 * Temporary response used for polygon hit detection.
 */
const T_RESPONSE = new Response();

/**
 * Tiny "point" polygon used for polygon hit detection.
 */
const TEST_POINT = new Box(new Vector(), 0.000001, 0.000001).toPolygon();

/// Helper Functions

/**
 * Flattens the specified array of points onto a unit vector axis,
 * resulting in a one dimensional range of the minimum and
 * maximum value on that axis.
 * @param points The points to flatten.
 * @param normal The unit vector axis to flatten on.
 * @param result An array.  After calling this function,
 *   result[0] will be the minimum value,
 *   result[1] will be the maximum value.
 */
function flattenPointsOn(points: Vector[], normal: Vector, result: number[]) {
    let min = Number.MAX_VALUE;
    let max = -Number.MAX_VALUE;
    const len = points.length;
    for (let i = 0; i < len; i++) {
        // The magnitude of the projection of the point onto the normal
        const dot = points[i].dot(normal);
        if (dot < min) { min = dot; }
        if (dot > max) { max = dot; }
    }
    result[0] = min;
    result[1] = max;
}

/**
 * Check whether two convex polygons are separated by the specified
 * axis (must be a unit vector).
 * @param aPos The position of the first polygon.
 * @param bPos The position of the second polygon.
 * @param aPoints The points in the first polygon.
 * @param bPoints The points in the second polygon.
 * @param axis The axis (unit sized) to test against.  The points of both polygons
 *   will be projected onto this axis.
 * @param response A Response object (optional) which will be populated
 *   if the axis is not a separating axis.
 * @returns true if it is a separating axis, false otherwise.  If false,
 *   and a response is passed in, information about how much overlap and
 *   the direction of the overlap will be populated.
 */
function isSeparatingAxis(aPos: Vector, bPos: Vector, aPoints: Vector[], bPoints: Vector[], axis: Vector, response?: Response): boolean {
    const rangeA = T_ARRAYS.pop()!;
    const rangeB = T_ARRAYS.pop()!;
    // The magnitude of the offset between the two polygons
    const offsetV = T_VECTORS.pop()!.copy(bPos).sub(aPos);
    const projectedOffset = offsetV.dot(axis);
    // Project the polygons onto the axis.
    flattenPointsOn(aPoints, axis, rangeA);
    flattenPointsOn(bPoints, axis, rangeB);
    // Move B's range to its position relative to A.
    rangeB[0] += projectedOffset;
    rangeB[1] += projectedOffset;
    // Check if there is a gap. If there is, this is a separating axis and we can stop
    if (rangeA[0] > rangeB[1] || rangeB[0] > rangeA[1]) {
        T_VECTORS.push(offsetV);
        T_ARRAYS.push(rangeA);
        T_ARRAYS.push(rangeB);
        return true;
    }
    // This is not a separating axis. If we're calculating a response, calculate the overlap.
    if (response) {
        let overlap = 0;
        // A starts further left than B
        if (rangeA[0] < rangeB[0]) {
            response.aInB = false;
            // A ends before B does. We have to pull A out of B
            if (rangeA[1] < rangeB[1]) {
                overlap = rangeA[1] - rangeB[0];
                response.bInA = false;
                // B is fully inside A.  Pick the shortest way out.
            } else {
                const option1 = rangeA[1] - rangeB[0];
                const option2 = rangeB[1] - rangeA[0];
                overlap = option1 < option2 ? option1 : -option2;
            }
            // B starts further left than A
        } else {
            response.bInA = false;
            // B ends before A ends. We have to push A out of B
            if (rangeA[1] > rangeB[1]) {
                overlap = rangeA[0] - rangeB[1];
                response.aInB = false;
                // A is fully inside B.  Pick the shortest way out.
            } else {
                const option1 = rangeA[1] - rangeB[0];
                const option2 = rangeB[1] - rangeA[0];
                overlap = option1 < option2 ? option1 : -option2;
            }
        }
        // If this is the smallest amount of overlap we've seen so far, set it as the minimum overlap.
        const absOverlap = Math.abs(overlap);
        if (absOverlap < response.overlap) {
            response.overlap = absOverlap;
            response.overlapN.copy(axis);
            if (overlap < 0) {
                response.overlapN.reverse();
            }
        }
    }
    T_VECTORS.push(offsetV);
    T_ARRAYS.push(rangeA);
    T_ARRAYS.push(rangeB);
    return false;
}

/**
 * Calculates which Voronoi region a point is on a line segment.
 * It is assumed that both the line and the point are relative to `(0,0)`
 *
 *            |       (0)      |
 *     (-1)  [S]--------------[E]  (1)
 *            |       (0)      |
 * @param line The line segment.
 * @param point The point.
 * @returns LEFT_VORONOI_REGION (-1) if it is the left region,
 *          MIDDLE_VORONOI_REGION (0) if it is the middle region,
 *          RIGHT_VORONOI_REGION (1) if it is the right region.
 */
function voronoiRegion(line: Vector, point: Vector): number {
    const len2 = line.len2();
    const dp = point.dot(line);
    // If the point is beyond the start of the line, it is in the
    // left voronoi region.
    if (dp < 0) { return LEFT_VORONOI_REGION; }
    // If the point is beyond the end of the line, it is in the
    // right voronoi region.
    else if (dp > len2) { return RIGHT_VORONOI_REGION; }
    // Otherwise, it's in the middle one.
    else { return MIDDLE_VORONOI_REGION; }
}

/// Constants for Voronoi regions

const LEFT_VORONOI_REGION = -1;
const MIDDLE_VORONOI_REGION = 0;
const RIGHT_VORONOI_REGION = 1;

/// Collision Tests

/**
 * Check if a point is inside a circle.
 * @param p The point to test.
 * @param c The circle to test.
 * @returns true if the point is inside the circle, false if it is not.
 */
function pointInCircle(p: Vector, c: Circle): boolean {
    const differenceV = T_VECTORS.pop()!.copy(p).sub(c.pos).sub(c.offset);
    const radiusSq = c.r * c.r;
    const distanceSq = differenceV.len2();
    T_VECTORS.push(differenceV);
    // If the distance between is smaller than the radius then the point is inside the circle.
    return distanceSq <= radiusSq;
}

/**
 * Check if a point is inside a convex polygon.
 * @param p The point to test.
 * @param poly The polygon to test.
 * @returns true if the point is inside the polygon, false if it is not.
 */
function pointInPolygon(p: Vector, poly: Polygon): boolean {
    TEST_POINT.pos.copy(p);
    T_RESPONSE.clear();
    let result = testPolygonPolygon(TEST_POINT, poly, T_RESPONSE);
    if (result) {
        result = T_RESPONSE.aInB;
    }
    return result;
}

/**
 * Check if two circles collide.
 * @param a The first circle.
 * @param b The second circle.
 * @param response Response object (optional) that will be populated if the circles intersect.
 * @returns true if the circles intersect, false if they don't.
 */
function testCircleCircle(a: Circle, b: Circle, response?: Response): boolean {
    // Check if the distance between the centers of the two
    // circles is greater than their combined radius.
    const differenceV = T_VECTORS.pop()!.copy(b.pos).add(b.offset).sub(a.pos).sub(a.offset);
    const totalRadius = a.r + b.r;
    const totalRadiusSq = totalRadius * totalRadius;
    const distanceSq = differenceV.len2();
    // If the distance is bigger than the combined radius, they don't intersect.
    if (distanceSq > totalRadiusSq) {
        T_VECTORS.push(differenceV);
        return false;
    }
    // They intersect.  If we're calculating a response, calculate the overlap.
    if (response) {
        const dist = Math.sqrt(distanceSq);
        response.a = a;
        response.b = b;
        response.overlap = totalRadius - dist;
        response.overlapN.copy(differenceV.normalize());
        response.overlapV.copy(differenceV).scale(response.overlap);
        response.aInB = a.r <= b.r && dist <= b.r - a.r;
        response.bInA = b.r <= a.r && dist <= a.r - b.r;
    }
    T_VECTORS.push(differenceV);
    return true;
}

/**
 * Check if a polygon and a circle collide.
 * @param polygon The polygon.
 * @param circle The circle.
 * @param response Response object (optional) that will be populated if they interset.
 * @returns true if they intersect, false if they don't.
 */
function testPolygonCircle(polygon: Polygon, circle: Circle, response?: Response): boolean {
    // Get the position of the circle relative to the polygon.
    const circlePos = T_VECTORS.pop()!.copy(circle.pos).add(circle.offset).sub(polygon.pos);
    const radius = circle.r;
    const radius2 = radius * radius;
    const points = polygon.calcPoints;
    const len = points.length;
    const edge = T_VECTORS.pop()!;
    const point = T_VECTORS.pop()!;

    // For each edge in the polygon:
    for (let i = 0; i < len; i++) {
        const next = i === len - 1 ? 0 : i + 1;
        const prev = i === 0 ? len - 1 : i - 1;
        let overlap = 0;
        let overlapN = null;

        // Get the edge.
        edge.copy(polygon.edges[i]);
        // Calculate the center of the circle relative to the starting point of the edge.
        point.copy(circlePos).sub(points[i]);

        // If the distance between the center of the circle and the point
        // is bigger than the radius, the polygon is definitely not fully in
        // the circle.
        if (response && point.len2() > radius2) {
            response.aInB = false;
        }

        // Calculate which Voronoi region the center of the circle is in.
        let region = voronoiRegion(edge, point);
        // If it's the left region:
        if (region === LEFT_VORONOI_REGION) {
            // We need to make sure we're in the RIGHT_VORONOI_REGION of the previous edge.
            edge.copy(polygon.edges[prev]);
            // Calculate the center of the circle relative the starting point of the previous edge
            const point2 = T_VECTORS.pop()!.copy(circlePos).sub(points[prev]);
            region = voronoiRegion(edge, point2);
            if (region === RIGHT_VORONOI_REGION) {
                // It's in the region we want.  Check if the circle intersects the point.
                const dist = point.len();
                if (dist > radius) {
                    // No intersection
                    T_VECTORS.push(circlePos);
                    T_VECTORS.push(edge);
                    T_VECTORS.push(point);
                    T_VECTORS.push(point2);
                    return false;
                } else if (response) {
                    // It intersects, calculate the overlap.
                    response.bInA = false;
                    overlapN = point.normalize();
                    overlap = radius - dist;
                }
            }
            T_VECTORS.push(point2);
            // If it's the right region:
        } else if (region === RIGHT_VORONOI_REGION) {
            // We need to make sure we're in the left region on the next edge
            edge.copy(polygon.edges[next]);
            // Calculate the center of the circle relative to the starting point of the next edge.
            point.copy(circlePos).sub(points[next]);
            region = voronoiRegion(edge, point);
            if (region === LEFT_VORONOI_REGION) {
                // It's in the region we want.  Check if the circle intersects the point.
                const dist = point.len();
                if (dist > radius) {
                    // No intersection
                    T_VECTORS.push(circlePos);
                    T_VECTORS.push(edge);
                    T_VECTORS.push(point);
                    return false;
                } else if (response) {
                    // It intersects, calculate the overlap.
                    response.bInA = false;
                    overlapN = point.normalize();
                    overlap = radius - dist;
                }
            }
            // Otherwise, it's the middle region:
        } else {
            // Need to check if the circle is intersecting the edge,
            // Change the edge into its "edge normal".
            const normal = edge.perp().normalize();
            // Find the perpendicular distance between the center of the
            // circle and the edge.
            const dist = point.dot(normal);
            const distAbs = Math.abs(dist);
            // If the circle is on the outside of the edge, there is no intersection.
            if (dist > 0 && distAbs > radius) {
                // No intersection
                T_VECTORS.push(circlePos);
                T_VECTORS.push(normal);
                T_VECTORS.push(point);
                return false;
            } else if (response) {
                // It intersects, calculate the overlap.
                overlapN = normal;
                overlap = radius - dist;
                // If the center of the circle is on the outside of the edge, or part of the
                // circle is on the outside, the circle is not fully inside the polygon.
                if (dist >= 0 || overlap < 2 * radius) {
                    response.bInA = false;
                }
            }
        }

        // If this is the smallest overlap we've seen, keep it.
        // (overlapN may be null if the circle was in the wrong Voronoi region).
        if (overlapN && response && Math.abs(overlap) < Math.abs(response.overlap)) {
            response.overlap = overlap;
            response.overlapN.copy(overlapN);
        }
    }

    // Calculate the final overlap vector - based on the smallest overlap.
    if (response) {
        response.a = polygon;
        response.b = circle;
        response.overlapV.copy(response.overlapN).scale(response.overlap);
    }
    T_VECTORS.push(circlePos);
    T_VECTORS.push(edge);
    T_VECTORS.push(point);
    return true;
}

/**
 * Check if a circle and a polygon collide.
 * **NOTE:** This is slightly less efficient than polygonCircle as it just
 * runs polygonCircle and reverses everything at the end.
 * @param circle The circle.
 * @param polygon The polygon.
 * @param response Response object (optional) that will be populated if they interset.
 * @returns true if they intersect, false if they don't.
 */
function testCirclePolygon(circle: Circle, polygon: Polygon, response?: Response): boolean {
    // Test the polygon against the circle.
    const result = testPolygonCircle(polygon, circle, response);
    if (result && response) {
        // Swap A and B in the response.
        const a = response.a;
        const aInB = response.aInB;
        response.overlapN.reverse();
        response.overlapV.reverse();
        response.a = response.b;
        response.b = a;
        response.aInB = response.bInA;
        response.bInA = aInB;
    }
    return result;
}

/**
 * Checks whether polygons collide.
 * @param a The first polygon.
 * @param b The second polygon.
 * @param response Response object (optional) that will be populated if they interset.
 * @returns true if they intersect, false if they don't.
 */
function testPolygonPolygon(a: Polygon, b: Polygon, response?: Response): boolean {
    const aPoints = a.calcPoints;
    const aLen = aPoints.length;
    const bPoints = b.calcPoints;
    const bLen = bPoints.length;
    let i;
    // If any of the edge normals of A is a separating axis, no intersection.
    for (i = 0; i < aLen; i++) {
        if (isSeparatingAxis(a.pos, b.pos, aPoints, bPoints, a.normals[i], response)) {
            return false;
        }
    }
    // If any of the edge normals of B is a separating axis, no intersection.
    for (i = 0; i < bLen; i++) {
        if (isSeparatingAxis(a.pos, b.pos, aPoints, bPoints, b.normals[i], response)) {
            return false;
        }
    }
    // Since none of the edge normals of A or B are a separating axis, there is an intersection
    // and we've already calculated the smallest overlap (in isSeparatingAxis).  Calculate the
    // final overlap vector.
    if (response) {
        response.a = a;
        response.b = b;
        response.overlapV.copy(response.overlapN).scale(response.overlap);
    }
    return true;
}

const SAT:Record<string, unknown> = {
    Vector,
    V: Vector,
    Box,
    Circle,
    Polygon,
    Response,
    isSeparatingAxis,
    pointInCircle,
    pointInPolygon,
    testCircleCircle,
    testPolygonCircle,
    testCirclePolygon,
    testPolygonPolygon
};
export default SAT;
