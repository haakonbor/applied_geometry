#ifndef CLOSED_SUBDIVISION_CURVE_H
#define CLOSED_SUBDIVISION_CURVE_H

#include <numeric>
#include <parametrics/gmpcurve.h>

namespace hbb {

template <typename T> class ClosedSubdivisionCurve : public GMlib::PCurve<T,3> {
    GM_SCENEOBJECT(ClosedSubdivisionCurve)

public:
    ClosedSubdivisionCurve(const GMlib::DVector<GMlib::Vector<T,3>>& p, int d);
    virtual ~ClosedSubdivisionCurve();

    //****************************************
    //****** Virtual public functions   ******
    //****************************************

    // from PCurve
    bool isClosed() const override;
    void sample(int m, int d = 0) override;

protected:
    // Virtual functions from PCurve, which have to be implemented locally
    void eval(T t, int d, bool l) const override;
    T getStartP() const override;
    T getEndP() const override;

    GMlib::DVector<GMlib::Vector<T, 3>> LaneRiesenfeldClosed(int k, int d);
    int  doublePart(GMlib::DVector<GMlib::Vector<T, 3>>& p, int n);
    void smoothPartClosed(GMlib::DVector<GMlib::Vector<T, 3>>& p, int n, int d);


    GMlib::DVector<GMlib::Vector<T, 3>> _p; // Points
    int _k = 3;     // Order
    int _d = 2;     // Degree

}; // END class ClosedSubdivisionCurve

//*****************************************
// Constructors and destructor           **
//*****************************************

template <typename T>
inline ClosedSubdivisionCurve<T>::ClosedSubdivisionCurve(const GMlib::DVector<GMlib::Vector<T,3>>& p, int d)
    : GMlib::PCurve<T, 3>(0,0,0), _p(p), _d(d), _k(d+1) {
}

template <typename T>
ClosedSubdivisionCurve<T>::~ClosedSubdivisionCurve() {}

//**************************************
//        Public local functons       **
//**************************************

//***************************************************
// Overrided (public) virtual functons from PCurve **
//***************************************************

template <typename T>
bool ClosedSubdivisionCurve<T>::isClosed() const { return true; }

template <typename T>
void ClosedSubdivisionCurve<T>::eval(T t, int d, bool /*l*/) const { }

//******************************************************
// Overrided (protected) virtual functons from PCurve **
//******************************************************

template <typename T>
T ClosedSubdivisionCurve<T>::getStartP() const { return T(0); }

template <typename T>
T ClosedSubdivisionCurve<T>::getEndP() const { return T(1); }


//**************************************
//        Protected local functons    **
//**************************************

template <typename T>
inline void ClosedSubdivisionCurve<T>::sample(int m, int d)
{
    // From original sample/presample functions
    this->_checkSampleVal(m, d);
    this->_visu.resize(1);

    // Algorithm 7 at page 112 in the book
    auto phi = LaneRiesenfeldClosed(m, _d);

    // Set the sample values to be the points found in the previous function
    _visu[0].sample_val.resize(phi.getDim());
    for (int i = 0; i < phi.getDim(); i++) {
        _visu[0].sample_val[i].setDim(1);
        _visu[0].sample_val[i][0] = phi[i];
    }

    // Fix the surrounding sphere
    _visu[0].sur_sphere.reset();
    computeSurroundingSphere(_visu[0].sample_val, _visu[0].sur_sphere);

    this->setEditDone();
}

template <typename T>
GMlib::DVector<GMlib::Vector<T, 3>> ClosedSubdivisionCurve<T>::LaneRiesenfeldClosed(int k, int d)
{
    int n = _p.getDim();                            // The number of intervals
    int m = (2 << k - 1) * n + 1;                   // The final number of points

    GMlib::DVector<GMlib::Vector<T, 3>> phi(m);     // The return vector - m points
    for (int i = 0; i < n; i++) {
        phi[i] = _p[i];                             // Inserting the initial points
    }

    phi[n] = _p[0];                                 // Closing the curve
    n++;

    for (int i = 0; i < k; i++) {                   // For each level of refinement
        n = doublePart(phi, n);
        smoothPartClosed(phi, n, d);
    }

    return phi;
}

template <typename T>
int ClosedSubdivisionCurve<T>::doublePart(GMlib::DVector<GMlib::Vector<T, 3>>& p, int n)
{
    for (int i = n - 1; i > 0; i--) {
        p[2 * i] = p[i];
        p[2 * i - 1] = (p[i] + p[i - 1]) / 2;
    }
    return 2 * n - 1;
}

template <typename T>
void ClosedSubdivisionCurve<T>::smoothPartClosed(GMlib::DVector<GMlib::Vector<T, 3>>& p, int n, int d)
{
    for (int j = 1; j < d; j++) {
        for (int i = 0; i < n - 1; i++) {
            p[i] = (p[i] + p[i + 1]) / 2;
        }
    }
    p[n - 1] = p[0];
}




}; // END namespace hbb

#endif // CLOSED_SUBDIVISION_CURVE_H
