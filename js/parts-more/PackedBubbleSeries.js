/**
 * (c) 2010-2018 Grzegorz Blachlinski, Sebastian Bochan
 *
 * License: www.highcharts.com/license
 */

'use strict';

import H from '../parts/Globals.js';
import '../parts/Utilities.js';
import '../parts/Axis.js';
import '../parts/Color.js';
import '../parts/Point.js';
import '../parts/Series.js';
import '../modules/networkgraph/networkgraph.src.js';


var seriesType = H.seriesType,
    defined = H.defined,
    pick = H.pick;


H.networkgraphIntegrations.packedbubble = {
    repulsiveForceFunction: function (d, k, node) {
        // todo: repNode
        // Used in API:
        // node.radius + repNode.radius?
        return Math.min(d, node.marker.radius + 5);
        // Force only for close nodes
    },
    barycenter: function () {
        var gravitationalConstant = this.options.gravitationalConstant,
            /* Add in next implementation
             *   xFactor = this.barycenter.xFactor,
             *   yFactor = this.barycenter.yFactor,
             */
            box = this.box,
            nodes = this.nodes;


        this.nodes.forEach(function (node) {
            if (!node.fixedPosition) {
                node.plotX -= (node.plotX - box.width / 2) *
                gravitationalConstant / (node.mass * Math.sqrt(nodes.length));
                node.plotY -= (node.plotY - box.height / 2) *
                gravitationalConstant / (node.mass * Math.sqrt(nodes.length));
            }
        });
    },
    repulsive: function (node, force, distanceXY, repNode) {
        var factor = force * this.diffTemperature / node.mass / node.degree;
        if (!node.fixedPosition) {
            node.plotX += distanceXY.x * factor;
            node.plotY += distanceXY.y * factor;
        }
        if (!repNode.fixedPosition) {
            repNode.plotX -= distanceXY.x * factor;
            repNode.plotY -= distanceXY.y * factor;
        }
    },
    integrate: function (layout, node) {
        /*
        Verlet without velocity:

            x(n+1) = 2 * x(n) - x(n-1) + A(T) * deltaT ^ 2

        where:
            - x(n+1) - new position
            - x(n) - current position
            - x(n-1) - previous position

        Assuming A(t) = 0 (no acceleration) and (deltaT = 1) we
        get:

            x(n+1) = x(n) + (x(n) - x(n-1))

        where:
            - (x(n) - x(n-1)) - position change

        TO DO:
        Consider Verlet with velocity to support additional
        forces. Or even Time-Corrected Verlet by Jonathan
        "lonesock" Dummer
        */
        var prevX = node.prevX,
            prevY = node.prevY,
            diffX = (node.plotX + node.dispX - prevX),
            diffY = (node.plotY + node.dispY - prevY);

        // Store for the next iteration:
        node.prevX = node.plotX + node.dispX;
        node.prevY = node.plotY + node.dispY;
        // Update positions, apply friction:
        node.plotX += diffX * -layout.options.friction;
        // /node.collisionNmb/node.neighbours;
        node.plotY += diffY * -layout.options.friction;
        // /node.collisionNmb/node.neighbours;

        node.temperature = layout.vectorLength({
            x: diffX * -layout.options.friction,
            y: diffY * -layout.options.friction
        });
    },
    getK: function (layout) {
        return Math.pow(
            layout.box.width * layout.box.height / layout.nodes.length,
            0.5
        );
    }
};

H.extend(
    H.layouts['reingold-fruchterman'].prototype,
    {

        setCircularPositions: function () {
            var box = this.box,
                nodes = this.nodes,
                nodesLength = nodes.length + 1,
                angle = 2 * Math.PI / nodesLength;


            // Initial positions are laid out along a small circle, appearing
            // as a cluster in the middle
            nodes.forEach(function (node, index) {
                node.plotX = node.prevX = pick(
                    node.plotX,
                    box.width / 2 + 80 * Math.cos(index * angle)
                );
                node.plotY = node.prevY = pick(
                    node.plotY,
                    box.height / 2 + 80 * Math.sin(index * angle)
                );

                node.dispX = 0;
                node.dispY = 0;
            });
        },
        getNeighbours: function () {
            var layout = this;
            layout.nodes.forEach(function (node) {
                node.neighbours = 1;
                layout.nodes.forEach(function (repNode) {
                    var distanceXY = layout.getDistXY(node, repNode);
                    var distanceR = layout.vectorLength(distanceXY);
                    if (
                        distanceR - (
                            repNode.marker.radius + node.marker.radius
                        ) < 1
                    ) {
                        node.neighbours += 0.01;
                    }
                });
            });
        },
        repulsiveForces: function () {
            var layout = this;
            layout.nodes.forEach(function (node) {
                node.degree = node.mass;
                layout.nodes.forEach(function (repNode) {
                    var force = 0,
                        distanceR,
                        distanceXY;

                    if (
                        // Node can not repulse itself:
                        node !== repNode &&
                        // Only close nodes affect each other:
                        /* layout.getDistR(node, repNode) < 2 * k && */
                        // Not dragged:
                        !node.fixedPosition
                    ) {
                        distanceXY = layout.getDistXY(node, repNode);
                        distanceR = (layout.vectorLength(distanceXY) -
                            (node.marker.radius + repNode.marker.radius + 5));

                        if (distanceR < 0) {
                            node.degree += 0.01;

                            // improvement needed v7.1.
                            if (layout.vectorLength(distanceXY) < 5) {
                                force = layout.repulsiveForce(
                                    -distanceR,
                                    layout.k,
                                    node,
                                    repNode
                                );
                            } else {
                                force = layout.repulsiveForce(
                                    -distanceR,
                                    layout.k,
                                    node,
                                    repNode
                                );
                            }
                        }

                        layout.force(
                            'repulsive',
                            node,
                            force * repNode.mass,
                            distanceXY,
                            repNode,
                            distanceR
                        );
                    }
                });
            });
        }
    }
);


/**
 * A packed bubble series is a two dimensional series type, where each point
 * renders a value in X, Y position. Each point is drawn as a bubble
 * where the bubbles don't overlap with each other and the radius
 * of the bubble related to the value.
 * Requires `highcharts-more.js`.
 *
 * @extends plotOptions.bubble
 * @excluding minSize,maxSize,connectNulls,keys,sizeByAbsoluteValue,
 * step,zMin,zMax,sizeBy,connectEnds
 * @product highcharts
 * @sample {highcharts} highcharts/demo/packed-bubble/
 *         Packed-bubble chart
 * @since 7.0.0
 * @optionparent plotOptions.packedbubble
 */

seriesType('packedbubble', 'bubble',
    {
        /**
         * Minimum bubble size. Bubbles will automatically size between the
         * `minSize` and `maxSize` to reflect the `z` value of each bubble.
         * Can be either pixels (when no unit is given), or a percentage of
         * the smallest one of the plot width and height.
         *
         * @type    {Number|String}
         * @sample  {highcharts} highcharts/plotoptions/bubble-size/ Bubble size
         * @since   3.0
         * @product highcharts highstock
         */
        minSize: '10%',
        /**
         * Maximum bubble size. Bubbles will automatically size between the
         * `minSize` and `maxSize` to reflect the `z` value of each bubble.
         * Can be either pixels (when no unit is given), or a percentage of
         * the smallest one of the plot width and height.
         *
         * @type    {Number|String}
         * @sample  {highcharts} highcharts/plotoptions/bubble-size/
         *          Bubble size
         * @since   3.0
         * @product highcharts highstock
         */
        maxSize: '50%',
        sizeBy: 'radius',
        zoneAxis: 'y',
        tooltip: {
            pointFormat: 'Value: {point.value}'
        },
        layoutAlgorithm: {
            /**
             * Repulsive force applied on a node. Passed are two arguments:
             * - `d` - which is current distance between two nodes
             * - `k` - which is desired distance between two nodes
             *
             * In `verlet` integration, defaults to:
             * `function (d, k) { return (k - d) / d * (k > d ? 1 : 0) }`
             *
             * @see         [layoutAlgorithm.integration](#series.networkgraph.layoutAlgorithm.integration)
             * @apioption   plotOptions.networkgraph.layoutAlgorithm.repulsiveForce
             * @sample      highcharts/series-networkgraph/forces/
             *              Custom forces with Euler integration
             * @sample      highcharts/series-networkgraph/cuboids/
             *              Custom forces with Verlet integration
             * @type        {Function}
             * @default function (d, k) { return k * k / d; }
             */

            /**
             * Attraction force applied on a node which is conected to another
             * node by a link. Passed are two arguments:
             * - `d` - which is current distance between two nodes
             * - `k` - which is desired distance between two nodes
             *
             * In `verlet` integration, defaults to:
             * `function (d, k) { return (k - d) / d; }`
             *
             * @see         [layoutAlgorithm.integration](#series.networkgraph.layoutAlgorithm.integration)
             * @apioption   plotOptions.networkgraph.layoutAlgorithm.attractiveForce
             * @sample      highcharts/series-networkgraph/forces/
             *              Custom forces with Euler integration
             * @sample      highcharts/series-networkgraph/cuboids/
             *              Custom forces with Verlet integration
             * @type        {Function}
             * @default function (d, k) { return k * k / d; }
             */

            /**
             * Ideal length (px) of the link between two nodes.
             * When not defined,
             * length is calculated as:
             * `Math.pow(availableWidth * availableHeight / nodesLength, 0.4);`
             *
             * Note: Because of the algorithm specification, length of each link
             * might be not exactly as specified.
             *
             * @type      {number}
             * @apioption series.networkgraph.layoutAlgorithm.linkLength
             * @sample    highcharts/series-networkgraph/styled-links/
             *            Numerical values
             * @defaults  undefined
             */

            /**
             * Initial layout algorithm for positioning nodes. Can be one of
             * built-in options ("circle", "random") or a function where
             * positions should be set on each node (`this.nodes`)
             *  as `node.plotX` and `node.plotY`
             *
             * @sample      highcharts/series-networkgraph/initial-positions/
             *              Initial positions with callback
             * @type        {String|Function}
             * @validvalue  ["circle", "random"]
             */
            initialPositions: 'circle',
            /**
             * Experimental. Enables live simulation of the algorithm
             * implementation. All nodes are animated as the forces applies on
             * them.
             *
             * @sample       highcharts/demo/network-graph/
             *               Live simulation enabled
             */
            enableSimulation: true,
            /**
             * Barnes-Hut approximation only.
             * Deteremines when distance between cell and node
             * is small enough to
             * caculate forces. Value of `theta` is compared
             * directly with quotient
             * `s / d`, where `s` is the size of the cell, and `d` is distance
             * between center of cell's mass and currently compared node.
             *
             * @see         [layoutAlgorithm.approximation](#series.networkgraph.layoutAlgorithm.approximation)
             * @since       7.1.0
             */
            theta: 0,
            /**
             * Approximation used to calculate repulsive forces affecting nodes.
             * By default, when calculateing net force, nodes are
             * compared against each other,
             * which gives O(N^2) complexity. Using Barnes-Hut
             * approximation, we decrease this to O(N log N), but the resulting
             * graph will have different layout.
             * Barnes-Hut approximation divides
             * space into rectangles via quad tree,
             * where forces exerted on nodes
             * are calculated directly for nearby cells, and for all others,
             * cells are treated as a separate node with center of mass.
             *
             * @see         [layoutAlgorithm.theta](#series.networkgraph.layoutAlgorithm.theta)
             * @validvalue  ["barnes-hut", "none"]
             * @since       7.1.0
             */
            approximation: 'none',
            /**
             * Type of the algorithm used when positioning nodes.
             *
             * @validvalue  ["reingold-fruchterman"]
             */
            type: 'reingold-fruchterman',
            /**
             * Integration type. Available options are `'euler'` and `'verlet'`.
             * Integration determines how forces are applied
             * on particles. In Euler
             * integration, force is applied direct as
             * `newPosition += velocity;`.
             * In Verlet integration, new position
             * is based on a previous posittion
             * without velocity:
             * `newPosition += previousPosition - newPosition`.
             *
             * Note that different integrations give different results as forces
             * are different.
             *
             * In Highcharts v7.0.x only `'euler'` integration was supported.
             *
             * @since       7.1.0
             * @sample      highcharts/series-networkgraph/forces/
             *              Custom forces with Euler integration
             * @validvalue  ["euler", "verlet"]
             */
            integration: 'packedbubble',
            /**
             * Max number of iterations before algorithm will stop. In general,
             * algorithm should find positions sooner, but when rendering huge
             * number of nodes, it is recommended to increase this value as
             * finding perfect graph positions can require more time.
             */
            maxIterations: 1000,
            /**
             * Gravitational const used in the barycenter
             * force of the algorithm.
             *
             * @sample      highcharts/series-networkgraph/forces/
             *              Custom forces with Euler integration
             */
            gravitationalConstant: 0.0325,
            /**
             * Friction applied on forces to prevent nodes rushing to fast to
             * the desired positions.
             */
            friction: -0.981
        }
    }, {
        forces: ['barycenter', 'repulsive'],
        pointArrayMap: ['value'],
        pointValKey: 'value',
        isCartesian: false,
        axisTypes: [],
        /**
         * Create a single array of all points from all series
         *
         * @param {Array} Array of all series objects
         * @return {Array} Returns the array of all points.
         *
         */
        accumulateAllPoints: function (series) {

            var chart = series.chart,
                allDataPoints = [],
                i, j;

            for (i = 0; i < chart.series.length; i++) {

                series = chart.series[i];

                if (series.visible || !chart.options.chart.ignoreHiddenSeries) {

                    // add data to array only if series is visible
                    for (j = 0; j < series.yData.length; j++) {
                        allDataPoints.push([
                            null, null,
                            series.yData[j],
                            series.index,
                            j,
                            {
                                id: j,
                                marker: {
                                    radius: 0
                                }
                            }
                        ]);
                    }
                }
            }

            return allDataPoints;
        },
        /* TODO
        init: function () {

            H.Series.prototype.init.apply(this, arguments);

            // When one series is modified, the others need to be recomputed
            H.addEvent(this, 'updatedData', function () {
                var self = this;
                this.chart.series.forEach(function (s) {
                    if (s.type === self.type) {
                        s.isDirty = true;
                    }
                });
            });

            return this;
        },
        */

        deferLayout: function () {
            var layoutOptions = this.options.layoutAlgorithm,
                points = this.points,
                graphLayoutsStorage = this.chart.graphLayoutsStorage,
                chartOptions = this.chart.options.chart,
                layout;

            this.nodes = points;
            this.nodes.forEach(function (node) {
                node.mass = Math.max(node.marker.radius / 10, 1);
                node.degree = 1;
                node.collisionNmb = 1;
            });

            if (!this.visible) {
                return;
            }

            if (!graphLayoutsStorage) {
                this.chart.graphLayoutsStorage = graphLayoutsStorage = {};
            }

            layout = graphLayoutsStorage[layoutOptions.type];

            if (!layout) {
                layoutOptions.enableSimulation =
                    !defined(chartOptions.forExport) ?
                        layoutOptions.enableSimulation :
                        !chartOptions.forExport;

                graphLayoutsStorage[layoutOptions.type] = layout =
                    new H.layouts[layoutOptions.type](layoutOptions);
            }

            this.layout = layout;

            layout.setArea(0, 0, this.chart.plotWidth, this.chart.plotHeight);
            layout.addSeries(this);
            layout.addNodes(this.nodes);
            layout.addLinks([]);
            this.points = points;
        },
        /**
         * Extend the base translate method to handle bubble size,
         * and correct positioning them
         */
        translate: function () {

            var series = this,
                chart = series.chart,
                data = series.data,
                index = series.index,
                point,
                radius,
                i;

            this.processedXData = this.xData;
            this.generatePoints();

            // merged data is an array with all of the data from all series
            if (!defined(chart.allDataPoints)) {
                chart.allDataPoints = series.accumulateAllPoints(series);
                // calculate radius for all added data
                series.getPointRadius();
            }

            // Set the shape type and arguments to be picked up in drawPoints
            for (i = 0; i < chart.allDataPoints.length; i++) {

                if (chart.allDataPoints[i][3] === index) {

                    // update the series points with the values from positions
                    // array
                    point = data[chart.allDataPoints[i][4]];
                    radius = chart.allDataPoints[i][2];
                    point.marker = H.extend(point.marker, {
                        radius: radius,
                        width: 2 * radius,
                        height: 2 * radius
                    });
                }
            }

            this.deferLayout();
        },
        /**
         * Check if two bubbles overlaps.
         * @param {Array} first bubble
         * @param {Array} second bubble
         *
         * @return {Boolean} overlap or not
         *
         */
        checkOverlap: function (bubble1, bubble2) {
            var diffX = bubble1[0] - bubble2[0], // diff of X center values
                diffY = bubble1[1] - bubble2[1], // diff of Y center values
                sumRad = bubble1[2] + bubble2[2]; // sum of bubble radius

            return (
                Math.sqrt(diffX * diffX + diffY * diffY) -
                Math.abs(sumRad)
            ) < -0.001;
        },
        /* Function that is adding one bubble based on positions and sizes
         * of two other bubbles, lastBubble is the last added bubble,
         * newOrigin is the bubble for positioning new bubbles.
         * nextBubble is the curently added bubble for which we are
         * calculating positions
         *
         * @param {Array} The closest last bubble
         * @param {Array} New bubble
         * @param {Array} The closest next bubble
         *
         * @return {Array} Bubble with correct positions
         *
         */
        positionBubble: function (lastBubble, newOrigin, nextBubble) {
            var sqrt = Math.sqrt,
                asin = Math.asin,
                acos = Math.acos,
                pow = Math.pow,
                abs = Math.abs,
                distance = sqrt( // dist between lastBubble and newOrigin
                    pow((lastBubble[0] - newOrigin[0]), 2) +
                    pow((lastBubble[1] - newOrigin[1]), 2)
                ),
                alfa = acos(
                    // from cosinus theorem: alfa is an angle used for
                    // calculating correct position
                    (
                        pow(distance, 2) +
                        pow(nextBubble[2] + newOrigin[2], 2) -
                        pow(nextBubble[2] + lastBubble[2], 2)
                    ) / (2 * (nextBubble[2] + newOrigin[2]) * distance)
                ),

                beta = asin( // from sinus theorem.
                    abs(lastBubble[0] - newOrigin[0]) /
                    distance
                ),
                // providing helping variables, related to angle between
                // lastBubble and newOrigin
                gamma = (lastBubble[1] - newOrigin[1]) < 0 ? 0 : Math.PI,
                // if new origin y is smaller than last bubble y value
                // (2 and 3 quarter),
                // add Math.PI to final angle

                delta = (lastBubble[0] - newOrigin[0]) *
                (lastBubble[1] - newOrigin[1]) < 0 ?
                    1 : -1, // (1st and 3rd quarter)
                finalAngle = gamma + alfa + beta * delta,
                cosA = Math.cos(finalAngle),
                sinA = Math.sin(finalAngle),
                posX = newOrigin[0] + (newOrigin[2] + nextBubble[2]) * sinA,
                // center of new origin + (radius1 + radius2) * sinus A
                posY = newOrigin[1] - (newOrigin[2] + nextBubble[2]) * cosA;

            return [
                posX,
                posY,
                nextBubble[2],
                nextBubble[3],
                nextBubble[4]
            ]; // the same as described before
        },
        /**
         * This is the main function responsible
         * for positioning all of the bubbles
         * allDataPoints - bubble array, in format [pixel x value,
         * pixel y value, radius,
         * related series index, related point index]
         *
         * @param {Array} All points from all series
         *
         * @return {Array} Positions of all bubbles
         *
         */
        placeBubbles: function (allDataPoints) {

            var series = this,
                checkOverlap = series.checkOverlap,
                positionBubble = series.positionBubble,
                bubblePos = [],
                stage = 1,
                j = 0,
                k = 0,
                calculatedBubble,
                sortedArr,
                arr = [],
                i;

            // sort all points
            sortedArr = allDataPoints.sort(function (a, b) {
                return b[2] - a[2];
            });

            if (sortedArr.length < 2) {
                // if length is 1,return only one bubble
                arr = [
                    0, 0,
                    sortedArr[0][0],
                    sortedArr[0][1],
                    sortedArr[0][2]
                ];
            } else {

                // create first bubble in the middle of the chart
                bubblePos.push([
                    [
                        0, // starting in 0,0 coordinates
                        0,
                        sortedArr[0][2], // radius
                        sortedArr[0][3], // series index
                        sortedArr[0][4]
                    ] // point index
                ]); // 0 level bubble

                bubblePos.push([
                    [
                        0,
                        0 - sortedArr[1][2] - sortedArr[0][2],
                        // move bubble above first one
                        sortedArr[1][2],
                        sortedArr[1][3],
                        sortedArr[1][4]
                    ]
                ]); // 1 level 1st bubble

                // first two already positioned so starting from 2
                for (i = 2; i < sortedArr.length; i++) {
                    sortedArr[i][2] = sortedArr[i][2] || 1;
                    // in case if radius is calculated as 0.

                    calculatedBubble = positionBubble(
                        bubblePos[stage][j],
                        bubblePos[stage - 1][k],
                        sortedArr[i]
                    ); // calculate initial bubble position

                    if (checkOverlap(calculatedBubble, bubblePos[stage][0])) {
                        /* if new bubble is overlapping with first bubble
                         * in current level (stage)
                         */

                        bubblePos.push([]);
                        k = 0;
                        /* reset index of bubble, used for
                         * positioning the bubbles
                         * around it, we are starting from first bubble in next
                         * stage because we are changing level to higher
                         */
                        bubblePos[stage + 1].push(
                            positionBubble(
                                bubblePos[stage][j],
                                bubblePos[stage][0],
                                sortedArr[i]
                            )
                        );
                        // (last added bubble, 1. from curr stage, new bubble)
                        stage++; // the new level is created, above current one
                        j = 0; // set the index of bubble in current level to 0
                    } else if (
                        stage > 1 && bubblePos[stage - 1][k + 1] &&
                        checkOverlap(
                            calculatedBubble, bubblePos[stage - 1][k + 1]
                        )
                    ) {
                        /* if new bubble is overlapping with one of the previous
                         * stage bubbles, it means that - bubble, used for
                         * positioning the bubbles around it has changed
                         * so we need to recalculate it
                         */
                        k++;
                        bubblePos[stage].push(
                            positionBubble(
                                bubblePos[stage][j],
                                bubblePos[stage - 1][k],
                                sortedArr[i]
                            )
                        );
                        // (last added bubble, prev stage bubble, new bubble)
                        j++;
                    } else { // simply add calculated bubble
                        j++;
                        bubblePos[stage].push(calculatedBubble);
                    }
                }
                series.chart.stages = bubblePos;
                // it may not be necessary but adding it just in case -
                // it is containing all of the bubble levels

                series.chart.rawPositions = [].concat.apply([], bubblePos);
                // bubble positions merged into one array
                arr = series.chart.rawPositions;
                series.resizeRadius();
            }
            return arr;
        },
        /**
         * The function responsible for resizing the bubble radius.
         * In shortcut: it is taking the initially
         * calculated positions of bubbles. Then it is calculating the min max
         * of both dimensions, creating something in shape of bBox.
         * The comparison of bBox and the size of plotArea
         * (later it may be also the size set by customer) is giving the
         * value how to recalculate the radius so it will match the size
         */
        resizeRadius: function () {

            var chart = this.chart,
                positions = chart.rawPositions,
                min = Math.min,
                max = Math.max,
                plotLeft = chart.plotLeft,
                plotTop = chart.plotTop,
                chartHeight = chart.plotHeight,
                chartWidth = chart.plotWidth,
                minX, maxX, minY, maxY,
                radius,
                bBox,
                spaceRatio,
                smallerDimension,
                i;

            minX = minY = Number.POSITIVE_INFINITY; // set initial values
            maxX = maxY = Number.NEGATIVE_INFINITY;

            for (i = 0; i < positions.length; i++) {
                radius = positions[i][2];
                minX = min(minX, positions[i][0] - radius);
                // (x center-radius) is the min x value used by specific bubble
                maxX = max(maxX, positions[i][0] + radius);
                minY = min(minY, positions[i][1] - radius);
                maxY = max(maxY, positions[i][1] + radius);
            }

            bBox = [maxX - minX, maxY - minY];
            spaceRatio = [
                (chartWidth - plotLeft) / bBox[0],
                (chartHeight - plotTop) / bBox[1]
            ];

            smallerDimension = min.apply([], spaceRatio);

            if (Math.abs(smallerDimension - 1) > 1e-10) {
                // if bBox is considered not the same width as possible size
                for (i = 0; i < positions.length; i++) {
                    positions[i][2] *= smallerDimension;
                }
                this.placeBubbles(positions);
            } else {
                /** if no radius recalculation is needed, we need to position
                 * the whole bubbles in center of chart plotarea
                 * for this, we are adding two parameters,
                 * diffY and diffX, that are related to differences
                 * between the initial center and the bounding box
                 */
                chart.diffY = chartHeight / 2 +
                    plotTop - minY - (maxY - minY) / 2;
                chart.diffX = chartWidth / 2 +
                    plotLeft - minX - (maxX - minX) / 2;
            }
        },

        /**
         * Calculate radius of bubbles in series.
         */
        getPointRadius: function () { // bubbles array

            var series = this,
                chart = series.chart,
                plotWidth = chart.plotWidth,
                plotHeight = chart.plotHeight,
                seriesOptions = series.options,
                smallestSize = Math.min(plotWidth, plotHeight),
                extremes = {},
                radii = [],
                allDataPoints = chart.allDataPoints,
                minSize,
                maxSize,
                value,
                radius;

            ['minSize', 'maxSize'].forEach(function (prop) {
                var length = parseInt(seriesOptions[prop], 10),
                    isPercent = /%$/.test(length);

                extremes[prop] = isPercent ?
                    smallestSize * length / 100 :
                    length;
            });

            chart.minRadius = minSize = extremes.minSize /
                Math.sqrt(allDataPoints.length);
            chart.maxRadius = maxSize = extremes.maxSize /
                Math.sqrt(allDataPoints.length);

            (allDataPoints || []).forEach(function (point, i) {
                value = point[2];

                radius = series.getRadius(
                    minSize,
                    maxSize,
                    minSize,
                    maxSize,
                    value
                );
                if (value === 0) {
                    radius = null;
                }
                allDataPoints[i][2] = radius;
                radii.push(radius);
            });

            this.radii = radii;
        },

        alignDataLabel: H.Series.prototype.alignDataLabel
    }, {
        getDegree: function () {
            return 1;
        }
    });

// Remove accumulated data points to redistribute all of them again
// (i.e after hiding series by legend)

H.addEvent(H.Chart, 'beforeRedraw', function () {
    if (this.allDataPoints) {
        delete this.allDataPoints;
    }
});

/**
 * A `packedbubble` series. If the [type](#series.packedbubble.type) option is
 * not specified, it is inherited from [chart.type](#chart.type).
 *
 * @type      {Object}
 * @extends   series,plotOptions.packedbubble
 * @excluding dataParser,dataURL,stack
 * @product   highcharts highstock
 * @apioption series.packedbubble
 */

/**
 * An array of data points for the series. For the `packedbubble` series type,
 * points can be given in the following ways:
 *
 * 1.  An array of `y` values.
 *
 *  ```js
 *     data: [5, 1, 20]
 *  ```
 *
 * 2.  An array of objects with named values. The objects are point
 * configuration objects as seen below. If the total number of data points
 * exceeds the series' [turboThreshold](#series.packedbubble.turboThreshold),
 * this option is not available.
 *
 *  ```js
 *     data: [{
 *         y: 1,
 *         name: "Point2",
 *         color: "#00FF00"
 *     }, {
 *         y: 5,
 *         name: "Point1",
 *         color: "#FF00FF"
 *     }]
 *  ```
 *
 * @type      {Array<Object|Array>}
 * @extends   series.line.data
 * @excluding marker
 * @sample    {highcharts} highcharts/series/data-array-of-objects/
 *            Config objects
 * @product   highcharts
 * @apioption series.packedbubble.data
 */

/**
 * @excluding enabled,enabledThreshold,height,radius,width
 * @apioption series.packedbubble.marker
 */
