Highcharts.chart('container', {
    chart: {
        type: 'packedbubble',
        height: '80%'
    },
    title: {
        text: 'TEAM TASKS'
    },
    tooltip: {
        useHTML: true,
        pointFormat: '<b>{point.name}:</b> {point.y}m CO<sub>2</sub>'
    },
    plotOptions: {
        packedbubble: {
            layoutAlgorithm: {
                mixSeries: false
            },
            dataLabels: {
                enabled: true,
                format: '{point.name}',
                filter: {
                    property: 'y',
                    operator: '>',
                    value: 250
                },
                style: {
                    color: 'black',
                    textOutline: 'none',
                    fontWeight: 'normal'
                }
            },
            minPointSize: 30
        }
    },
    series: [{
        link: {
            width: 0
        },
        name: 'Josephine',
        data: (function () {
            var d = [],
                points = 25;
            for (var i = 0; i < points; i++) {
                d.push([i, i * 10]);
            }
            return d;
        }())
    }, {
        name: 'Sophie',
        link: {
            width: 0
        },
        data: (function () {
            var d = [],
                points = 50;
            for (var i = 0; i < points; i++) {
                d.push([i, i * 10]);
            }
            return d;
        }())
    }, {
        name: 'Michael',
        link: {
            width: 0
        },
        data: (function () {
            var d = [],
                points = 20;
            for (var i = 0; i < points; i++) {
                d.push([i, i * 10]);
            }
            return d;
        }())
    }, {
        name: 'Jon',
        link: {
            width: 0
        },
        data: (function () {
            var d = [],
                points = 20;
            for (var i = 0; i < points; i++) {
                d.push([i, i * 10]);
            }
            return d;
        }())
    }]
});
