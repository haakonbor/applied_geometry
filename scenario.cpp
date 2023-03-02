
#include <iostream>

#include "application/b_spline_curve.h"
#include "application/blending_spline_curve.h"
#include "application/blending_spline_surface.h"
#include "application/closed_subdivision_curve.h"
#include "application/model_curve1.h"
#include "parametrics/curves/gmpcircle.h"
#include "parametrics/surfaces/gmpcylinder.h"
#include "parametrics/surfaces/gmpplane.h"
#include "parametrics/surfaces/gmptorus.h"
#include "scenario.h"
#include "testtorus.h"

// hidmanager
#include "hidmanager/defaulthidmanager.h"

// gmlib
#include <scene/light/gmpointlight.h>
#include <scene/sceneobjects/gmpathtrack.h>
#include <scene/sceneobjects/gmpathtrackarrows.h>

// qt
#include <QQuickItem>

template <typename T>
inline std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  out << v.size() << std::endl;
  for (uint i = 0; i < v.size(); i++)
    out << " " << v[i];
  out << std::endl;
  return out;
}

void Scenario::initializeScenario() {

  // Insert a light
  GMlib::Point<GLfloat, 3> init_light_pos(2.0, 4.0, 10);
  GMlib::PointLight *light =
      new GMlib::PointLight(GMlib::GMcolor::white(), GMlib::GMcolor::white(),
                            GMlib::GMcolor::white(), init_light_pos);
  light->setAttenuation(0.8f, 0.002f, 0.0008f);
  this->scene()->insertLight(light, false);

  // Insert Sun
  this->scene()->insertSun();

  // Default camera parameters
  int init_viewport_size = 600;
  GMlib::Point<float, 3> init_cam_pos(0.0f, 0.0f, 0.0f);
  GMlib::Vector<float, 3> init_cam_dir(0.0f, 1.0f, 0.0f);
  GMlib::Vector<float, 3> init_cam_up(1.0f, 0.0f, 0.0f);

  // Projection cam
  auto proj_rcpair = createRCPair("Projection");
  proj_rcpair.camera->set(init_cam_pos, init_cam_dir, init_cam_up);
  proj_rcpair.camera->setCuttingPlanes(1.0f, 8000.0f);
  proj_rcpair.camera->rotateGlobal(GMlib::Angle(-45),
                                   GMlib::Vector<float, 3>(1.0f, 0.0f, 0.0f));
  proj_rcpair.camera->translateGlobal(
      GMlib::Vector<float, 3>(0.0f, -20.0f, 20.0f));
  scene()->insertCamera(proj_rcpair.camera.get());
  proj_rcpair.renderer->reshape(
      GMlib::Vector<int, 2>(init_viewport_size, init_viewport_size));

  /***************************************************************************
   *                                                                         *
   * Standar example, including path track and path track arrows             *
   *                                                                         *
   ***************************************************************************/

  GMlib::Material mm(GMlib::GMmaterial::polishedBronze());
  mm.set(45.0);

  GMlib::DVector<GMlib::Vector<float, 3>> c(8);
  //  c[0] = {0, -5, 0};
  //  c[1] = {-5, -5, -2};
  //  c[2] = {-5, 0, -4};
  //  c[3] = {-5, 5, -6};
  //  c[4] = {0, 5, -8};
  //  c[5] = {5, 5, -10};
  //  c[6] = {5, 0, -12};
  //  c[7] = {5, -5, -10};
  //  c[8] = {0, -5, -8};
  //  c[9] = {-5, -5, -6};
  //  c[10] = {-5, 0, -4};
  //  c[11] = {-5, 5, -2};
  //  c[12] = {0, 0, 0};

  c[0] = {0, 0, 0};
  c[1] = {5, 0, 0};
  c[2] = {10, 5, 0};
  c[3] = {10, 10, 0};
  c[4] = {5, 15, 0};
  c[5] = {10, 20, 0};
  c[6] = {5, 25, 0};
  c[7] = {0, 20, 0};

  /* DEFAULT DEMO */
  //  auto ptom = new TestTorus(1.0f, 0.4f, 0.6f);
  //  ptom->toggleDefaultVisualizer();
  //  ptom->sample(60,60,1,1);
  //  this->scene()->insert(ptom);
  //  auto ptrack = new GMlib::PathTrack();
  //  ptrack->setLineWidth(2);
  //  ptom->insert(ptrack);
  //  auto ptrack2 = new GMlib::PathTrackArrows();
  //  ptrack2->setArrowLength(2);
  //  ptom->insert(ptrack2);

  /* 2ND DEGREE B-SPLINE */
  //  auto bspline = new hbb::My_B_spline<float>(c);
  //  bspline->toggleDefaultVisualizer();
  //  bspline->sample(60, 0);
  //  bspline->setLineWidth(4);
  //  this->scene()->insert(bspline);

  /* CLOSED SUBDIVISION CURVE */
  //  // 2nd-degree curve
  //  auto closed_subdivision_curve1 = new hbb::ClosedSubdivisionCurve<float>(c,
  //  2); closed_subdivision_curve1->toggleDefaultVisualizer();
  //  closed_subdivision_curve1->sample(3, 0);
  //  closed_subdivision_curve1->setLineWidth(4);
  //  closed_subdivision_curve1->setColor(GMlib::GMcolor::red());
  //  this->scene()->insert(closed_subdivision_curve1);

  //  // 3rd-degree  curve
  //  auto closed_subdivision_curve2 = new hbb::ClosedSubdivisionCurve<float>(c,
  //  3); closed_subdivision_curve2->toggleDefaultVisualizer();
  //  closed_subdivision_curve2->sample(3, 0);
  //  closed_subdivision_curve2->setLineWidth(4);
  //  closed_subdivision_curve2->setColor(GMlib::GMcolor::blue());
  //  this->scene()->insert(closed_subdivision_curve2);

  /* BLENDING SPLINE CURVE */
  /* MODEL CURVE */
  auto curve = new hbb::Model_curve1<float>();
  curve->toggleDefaultVisualizer();
  curve->sample(100, 1);
  curve->setLineWidth(1);
  this->scene()->insert(curve);
  auto blending_spline_curve = new hbb::Blending_spline_curve<float>(curve, 4);
  blending_spline_curve->toggleDefaultVisualizer();
  blending_spline_curve->sample(100, 0);
  blending_spline_curve->setLineWidth(4);
  blending_spline_curve->setColor(GMlib::GMcolor::blue());
  this->scene()->insert(blending_spline_curve);

  /* BLENDING SPLINE SURFACE */
  //  auto cylinder = new GMlib::PCylinder<float>();
  //  auto torus = new GMlib::PTorus<float>();
  //  auto plane =
  //      new GMlib::PPlane<float>(GMlib::Point<float, 3>(0.0f, 0.0f, 0.0f),
  //                               GMlib::Point<float, 3>(8.0f, 0.0f, 0.0f),
  //                               GMlib::Point<float, 3>(0.0f, 8.0f, 0.0f));

  //  // choose between cylinder, torus and plane
  //  auto model_surface = torus;

  //  auto blending_spline_surface =
  //      new hbb::Blending_spline_surface<float>(model_surface, 4, 4);
  //  blending_spline_surface->toggleDefaultVisualizer();
  //  blending_spline_surface->sample(50, 50, 1, 1);
  //  this->scene()->insert(blending_spline_surface);
}

void Scenario::cleanupScenario() {}

void Scenario::callDefferedGL() {

  GMlib::Array<const GMlib::SceneObject *> e_obj;
  this->scene()->getEditedObjects(e_obj);

  for (int i = 0; i < e_obj.getSize(); i++)
    if (e_obj(i)->isVisible())
      e_obj[i]->replot();
}
