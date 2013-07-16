
// From GLE
// http://developer.apple.com/library/mac/#samplecode/glut/Listings/gle_vvector_h.html

/* ========================================================== */
/* inverse of matrix 
 *
 * Compute inverse of matrix a, returning determinant m and 
 * inverse b
 */
 
#define INVERT_3X3(b,det,a)         \
    {                       \
    double tmp;                  \
    DETERMINANT_3X3 (det, a);            \
    tmp = 1.0 / (det);               \
    SCALE_ADJOINT_3X3 (b, tmp, a);       \
    }
