const symbolExx = {
  draw: function(context, size) {
    var n = 5;
    var m = n - 1;
    var r = Math.sqrt(size / Math.PI)/n; // to match radius of symbolCircle
    context.moveTo(-m * r, n * r);
    context.lineTo(0, r);
    context.lineTo(m * r, n * r);
    context.lineTo(n * r, m * r);
    context.lineTo(r, 0);
    context.lineTo(n * r, -m * r);
    context.lineTo(m * r, -n * r);
    context.lineTo(0, -r);
    context.lineTo(-m * r, -n * r);
    context.lineTo(-n * r, -m * r);
    context.lineTo(-r, 0);
    context.lineTo(-n * r, m * r);
    context.closePath();
  }
};
