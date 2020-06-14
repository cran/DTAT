const protoTHplot = { 
  width: 0.25 * width, 
  height: height * 0.45, 
  margin: dsMargin,
  ndoses: data.doses.length,
};

function thFactory(opts, proto = protoTHplot) {
  const thPlot = Object.assign({}, proto, opts);

  const margin = thPlot.margin;
  margin.top = 0;
  const horizMargins = margin.left + margin.right;
  const vertMargins = margin.top + margin.bottom;
  thPlot.svg = div.select('#ds-plot')
      .append('svg')
      .attr('id', thPlot.id || 'thPlot')
      .attr('viewBox',`0 0 ${thPlot.width+horizMargins} ${thPlot.height+vertMargins}`)
      .attr('width', thPlot.width + horizMargins) 
      .attr('height', thPlot.height + vertMargins)
    .append('g')
      //.attr('id','container')
      .attr('transform', `translate(${margin.left}, ${margin.top})`);

  return thPlot;
}

function renderTHplot(opts) {
  const thPlot = thFactory(opts);

  // Construct a curvilinear scale using the dose-survival data
  const dose_surv_map = data[data.length - 1].map(x => x.surv).slice(0,-1);
  dose_surv_map.unshift(1.0); // map 'dose 0' to '100% dose-survival'
  const x = d3.scaleLinear()
    .domain(d3.range(dose_surv_map.length))
    .range(dose_surv_map.map(f => thPlot.width*(1-f)));

  const dsteps = thPlot.ndoses - 1;

  const y = d3.scaleLinear()
    .range([0.9*thPlot.height, 0.1*thPlot.height])
    .domain([0, 3]);
  
  const xAxis = d3.axisBottom().scale(x);
  xAxis.tickValues(d3.range(1, dose_surv_map.length));
  xAxis.tickFormat(d3.format('.0f'));
  const yAxis = d3.axisRight().scale(y);
  yAxis.tickValues(d3.range(0,thPlot.ndoses+1));
  yAxis.tickFormat(d3.format('.0f'));

  thPlot.svg.append('g')
    .attr('class','axis')
    .attr('transform', `translate(0, ${thPlot.height})`)
    .call(xAxis);

  thPlot.svg.append('g')
    .attr('class','axis')
    .attr('transform', `translate(${thPlot.width}, 0)`)
    .call(yAxis);

  thPlot.svg.append("text")
    .attr("class", "x label")
    .attr("text-anchor", "middle")
    .attr("x", thPlot.width / 2)
    .attr("y", thPlot.height + 30)
    .text("Dose Level");

  thPlot.svg.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "middle")
    .attr("x", -thPlot.height*0.5)
    .attr("y", thPlot.width + 25)
    .attr("dy", ".75em")
    .attr("transform", "rotate(-90)")
    .text("On-target Toxicity Score");

  swimPlot.responses.forEach(i => {
    thPlot.svg.append('path') // TODO: vary symbol by response
        .datum({id: i.id, dose: i.dose, response: i.response, ott: i.ott})
        .attr("class", "response")
        .attr('d', responseMarker(i).size(thPlot.width*0.2))
        .attr('stroke', "black")
        .attr('transform', `translate(${x(i.dose)}, ${y(i.ott)})`)
        .attr("assessment", i.response);
  });

  return thPlot;
}
