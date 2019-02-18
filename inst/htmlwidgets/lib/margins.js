// Size each of the 2 DIVs appropriately...
const makeMargin = function(left, right) {
    return {top: 30,
            right: right,
            bottom: 50,
            left: left};
};

const oxMargin = makeMargin(40, 60);
const dsMargin = makeMargin(20, 50);

const horiz = m => m.left + m.right;
const vert = m => m.top + m.bottom;

width = width - horiz(oxMargin) - horiz(dsMargin);
height = height - vert(oxMargin);
