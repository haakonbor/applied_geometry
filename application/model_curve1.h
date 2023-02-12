#ifndef MODEL_CURVE1_H
#define MODEL_CURVE1_H

#include <parametrics/gmpcurve.h>

namespace hbb {

template <typename T> class Model_curve1 : public GMlib::PCurve<T, 3> {
  GM_SCENEOBJECT(Model_curve1)

public:
  Model_curve1();
  Model_curve1(const Model_curve1<T> &copy);
  virtual ~Model_curve1();

  //****************************************
  //****** Virtual public functions   ******
  //****************************************

  // from PCurve
  bool isClosed() const override;

protected:
  // Virtual functions from PCurve, which have to be implemented locally
  void eval(T t, int d, bool l) const override;
  T getStartP() const override;
  T getEndP() const override;

}; // END class Model_curve1

//*****************************************
// Constructors and destructor           **
//*****************************************

/*! Model_curve1<T>::Model_curve1( )
 *  Default constructor to make a curve.
 *
 */
template <typename T>
inline Model_curve1<T>::Model_curve1() : GMlib::PCurve<T, 3>(20, 0, 1) {}

/*! Model_curve1<T>::Model_curve1(const Model_curve1<T>& copy )
 *  A copy constructor
 *  Making a copy of the curve
 *
 *  \param[in] copy The curve to copy
 */
template <typename T>
inline Model_curve1<T>::Model_curve1(const Model_curve1<T> &copy)
    : GMlib::PCurve<T, 3>(copy) {}

/*! Model_curve1<T>::~Model_curve1()
 *  The destructor
 *  clean up and destroy all private data
 */
template <typename T> Model_curve1<T>::~Model_curve1() {}

//**************************************
//        Public local functons       **
//**************************************

//***************************************************
// Overrided (public) virtual functons from PCurve **
//***************************************************

/*! bool Model_curve1<T>::isClosed() const
 *  To tell that this curve (circle) is closed.
 */
template <typename T> bool Model_curve1<T>::isClosed() const { return false; }

//******************************************************
// Overrided (protected) virtual functons from PCurve **
//******************************************************

/*! void Model_curve1<T>::eval( T t, int d, bool l ) const
 *  Evaluation of the curve at a given parameter value
 *  To compute position and d derivatives at parameter value t on the curve.
 *  7 derivatives are implemented
 *
 *  \param  t[in]  The parameter value to evaluate at
 *  \param  d[in]  The number of derivatives to compute
 *  \param  l[in]  (dummy) because left and right are always equal
 */
template <typename T> void Model_curve1<T>::eval(T t, int d, bool /*l*/) const {

  this->_p.setDim(d + 1);

  this->_p[0][0] = std::cos(3 * t);
  this->_p[0][1] = std::sin(t);
  this->_p[0][2] = std::cos(2 * t);

  if (this->_dm == GMlib::GM_DERIVATION_EXPLICIT) {
    if (d > 0) {
      this->_p[1][0] = 3 * std::sin(3 * t);
      this->_p[1][1] = -std::cos(t);
      this->_p[1][2] = 2 * std::sin(2 * t);
    }
  }
}

/*! T Model_curve1<T>::getStartP() const
 *  Provides the start parameter value associated with
 *  the eval() function implemented above.
 *  (the start parameter value = 0).
 *
 *  \return The parametervalue at start of the internal domain
 */
template <typename T> T Model_curve1<T>::getStartP() const { return T(0); }

/*! T Model_curve1<T>::getEndP() const
 *  Provides the end parameter value associated with
 *  the eval() function implemented above.
 *  (the end parameter value = 2*Pi).
 *
 *  \return The parametervalue at end of the internal domain
 */
template <typename T> T Model_curve1<T>::getEndP() const { return T(M_2PI); }

} // namespace hbb

#endif // MODEL_CURVE1_H
