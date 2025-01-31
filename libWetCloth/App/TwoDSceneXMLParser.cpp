//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and
// Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TwoDSceneXMLParser.h"

#include <fstream>
#include <string>

#include "MathDefs.h"

void TwoDSceneXMLParser::loadExecutableSimulation(
    const std::string& file_name, bool rendering_enabled,
    std::shared_ptr<ParticleSimulation>& execsim, Camera& cam, scalar& dt,
    scalar& max_time, scalar& steps_per_sec_cap, renderingutils::Color& bgcolor,
    std::string& description, std::string& scenetag, bool& cam_inited,
    const std::string& input_bin) {
  // Load the xml document
  std::vector<char> xmlchars;
  rapidxml::xml_document<> doc;
  loadXMLFile(file_name, xmlchars, doc);

  // Attempt to locate the root node
  rapidxml::xml_node<>* node = doc.first_node("scene");
  if (node == NULL) {
    std::cerr << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse xml scene file. Failed to locate root "
                 "<scene> node. Exiting."
              << std::endl;
    exit(1);
  }

  // Determine what simulation type this is (particle, rigid body, etc)
  std::string simtype;
  loadSimulationType(node, simtype);

  // Parse common state
  loadMaxTime(node, max_time);
  loadMaxSimFrequency(node, steps_per_sec_cap);
  loadBackgroundColor(node, bgcolor);
  loadSceneDescriptionString(node, description);
  loadSceneTag(node, scenetag);

  cam_inited = loadCamera(node, cam);

  // Parse the user-requested simulation type. The default is a particle
  // simulation.
  loadParticleSimulation(rendering_enabled, execsim, dt, bgcolor, node,
                         input_bin);
}

void TwoDSceneXMLParser::loadBucketInfo(
    rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene) {
  rapidxml::xml_node<>* nd = node->first_node("bucketinfo");
  scalar bucket_size = 5.0;
  int num_cells = 5;
  int kernel_order = 2;

  if (nd) {
    if (nd->first_attribute("size")) {
      std::string attribute(nd->first_attribute("size")->value());
      if (!stringutils::extractFromString(attribute, bucket_size)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of size for bucketinfo. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("numcells")) {
      std::string attribute(nd->first_attribute("numcells")->value());
      if (!stringutils::extractFromString(attribute, num_cells)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of numcells for bucketinfo. Value "
                     "must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("kernelorder")) {
      std::string attribute(nd->first_attribute("kernelorder")->value());
      if (!stringutils::extractFromString(attribute, kernel_order)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of numcells for bucketinfo. Value "
                     "must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }
  }

  twodscene->setBucketInfo(bucket_size, num_cells, kernel_order);
}

void TwoDSceneXMLParser::loadParticleSimulation(
    bool rendering_enabled, std::shared_ptr<ParticleSimulation>& execsim,
    scalar& dt, renderingutils::Color& bgcolor, rapidxml::xml_node<>* node,
    const std::string& input_bin) {
  auto scene = std::make_shared<TwoDScene>();

  // Integrator/solver
  std::shared_ptr<SceneStepper> scene_stepper = NULL;
  loadIntegrator(node, scene_stepper, dt);
  assert(scene_stepper != NULL);
  assert(dt > 0.0);

  // Scene
  loadSimInfo(node, scene);
  loadBucketInfo(node, scene);

  int mg_part = 0, mg_df = 0;
  loadParticles(node, scene, mg_part);

  loadElasticParameters(node, scene, dt);

  int maxgroup = std::max(mg_part, mg_df);
  scene->resizeGroups(0);
  scene->updateRestPos();
  scene->initGroupPos();

  loadHairs(node, scene, dt);
  loadHairPose(node, scene);
  loadScripts(node, scene);

  scene->initGaussSystem();
  scene->updateShapeFactor();
  scene->updateParticleBoundingBox();
  scene->rebucketizeParticles();
  scene->resampleNodes();
  scene->updateManifoldOperators();
  scene->computeWeights(0.0);
  scene->updatePlasticity(0.0);
  scene->computedEdFe();
  scene->updateOrientation();

  scene->updateSolidPhi();
  scene->updateSolidWeights();
  scene->updateIntersection();
  scene->mapParticleNodesAPIC();

  scene->saveParticleVelocity();

  // Forces
  scene->loadAttachForces();
  loadSpringForces(node, scene);
  loadSimpleGravityForces(node, scene);

  std::shared_ptr<TwoDSceneRenderer> scene_renderer = NULL;
  if (rendering_enabled) {
    scene_renderer = std::make_shared<TwoDSceneRenderer>(*scene);
    scene_renderer->updateParticleSimulationState(*scene);
  }

  execsim = std::make_shared<ParticleSimulation>(scene, scene_stepper,
                                                 scene_renderer);

  if (!input_bin.empty()) {
    execsim->readPos(input_bin);
    scene->updateGaussSystem(0.0);
  }
}

void TwoDSceneXMLParser::loadSpringForces(
    rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene) {
  assert(node != NULL);

  int forcenum = 0;
  // Extract the edge the force acts across
  for (rapidxml::xml_node<>* nd = node->first_node("springforce"); nd;
       nd = nd->next_sibling("springforce")) {
    int edge = -1;

    if (nd->first_attribute("edge")) {
      std::string attribute(nd->first_attribute("edge")->value());
      if (!stringutils::extractFromString(attribute, edge)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of edge attribute for springforce "
                  << forcenum << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of edge attribute for springforce "
                << forcenum << ". Exiting." << std::endl;
      exit(1);
    }

    // Extract the spring stiffness
    scalar k = std::numeric_limits<scalar>::signaling_NaN();
    if (nd->first_attribute("k")) {
      std::string attribute(nd->first_attribute("k")->value());
      if (!stringutils::extractFromString(attribute, k)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of k attribute for springforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse k attribute for springforce " << forcenum
                << ". Exiting." << std::endl;
      exit(1);
    }

    // Extract the spring rest length
    scalar l0 = std::numeric_limits<scalar>::signaling_NaN();
    if (nd->first_attribute("l0")) {
      std::string attribute(nd->first_attribute("l0")->value());
      if (!stringutils::extractFromString(attribute, l0)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of l0 attribute for springforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse l0 attribute for springforce " << forcenum
                << ". Exiting." << std::endl;
      exit(1);
    }

    // Extract the optional damping coefficient
    scalar b = 0.0;
    if (nd->first_attribute("b")) {
      std::string attribute(nd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attribute, b)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of b attribute for springforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    // std::cout << "Springforce: " << forcenum << "    i: " << newedge.first <<
    // "   j: " << newedge.second << "   k: " << k << "   l0: " << l0 <<
    // std::endl;

    twodscene->insertForce(
        std::make_shared<SpringForce>(twodscene->getEdge(edge), k, l0, b));

    forcenum++;
  }

  // SpringForce( const std::pair<int,int>& endpoints, const scalar& k, const
  // scalar& l0 )
}

void TwoDSceneXMLParser::loadSimpleGravityForces(
    rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene) {
  assert(node != NULL);

  // Load each constant force
  int forcenum = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("simplegravity"); nd;
       nd = nd->next_sibling("simplegravity")) {
    Vector3s constforce;
    constforce.setConstant(0.0);

    // Extract the x component of the force
    if (nd->first_attribute("fx")) {
      std::string attribute(nd->first_attribute("fx")->value());
      if (!stringutils::extractFromString(attribute, constforce.x())) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fx attribute for constantforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    // Extract the y component of the force
    if (nd->first_attribute("fy")) {
      std::string attribute(nd->first_attribute("fy")->value());
      if (!stringutils::extractFromString(attribute, constforce.y())) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fy attribute for constantforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    // Extract the z component of the force
    if (nd->first_attribute("fz")) {
      std::string attribute(nd->first_attribute("fz")->value());
      if (!stringutils::extractFromString(attribute, constforce.z())) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fz attribute for constantforce "
                  << forcenum << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    twodscene->insertForce(std::make_shared<SimpleGravityForce>(constforce));

    ++forcenum;
  }
}

bool TwoDSceneXMLParser::loadCamera(rapidxml::xml_node<>* node,
                                    Camera& camera) {
  rapidxml::xml_node<>* nd = node->first_node("camera");
  if (nd) {
    Eigen::Quaterniond rotation(1, 0, 0, 0);
    rapidxml::xml_node<>* nd_rot = nd->first_node("rotation");

    if (nd_rot) {
      if (nd_rot->first_attribute("x")) {
        std::string attribute(nd_rot->first_attribute("x")->value());
        if (!stringutils::extractFromString(attribute, rotation.x())) {
          std::cerr << outputmod::startred
                    << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                    << " Failed to parse value of x attribute for rotation. "
                       "Value must be scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_rot->first_attribute("y")) {
        std::string attribute(nd_rot->first_attribute("y")->value());
        if (!stringutils::extractFromString(attribute, rotation.y())) {
          std::cerr << outputmod::startred
                    << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                    << " Failed to parse value of y attribute for rotation. "
                       "Value must be scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_rot->first_attribute("z")) {
        std::string attribute(nd_rot->first_attribute("z")->value());
        if (!stringutils::extractFromString(attribute, rotation.z())) {
          std::cerr << outputmod::startred
                    << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                    << " Failed to parse value of z attribute for rotation. "
                       "Value must be scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_rot->first_attribute("w")) {
        std::string attribute(nd_rot->first_attribute("w")->value());
        if (!stringutils::extractFromString(attribute, rotation.w())) {
          std::cerr << outputmod::startred
                    << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                    << " Failed to parse value of w attribute for rotation. "
                       "Value must be scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
    }

    Eigen::Vector3d center(0, 0, 0);
    rapidxml::xml_node<>* nd_center = nd->first_node("center");

    if (nd_center) {
      if (nd_center->first_attribute("x")) {
        std::string attribute(nd_center->first_attribute("x")->value());
        if (!stringutils::extractFromString(attribute, center.x())) {
          std::cerr << outputmod::startred
                    << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                    << " Failed to parse value of x attribute for center. "
                       "Value must be scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_center->first_attribute("y")) {
        std::string attribute(nd_center->first_attribute("y")->value());
        if (!stringutils::extractFromString(attribute, center.y())) {
          std::cerr << outputmod::startred
                    << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                    << " Failed to parse value of y attribute for center. "
                       "Value must be scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
      if (nd_center->first_attribute("z")) {
        std::string attribute(nd_center->first_attribute("z")->value());
        if (!stringutils::extractFromString(attribute, center.z())) {
          std::cerr << outputmod::startred
                    << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                    << " Failed to parse value of z attribute for center. "
                       "Value must be scalar. Exiting."
                    << std::endl;
          exit(1);
        }
      }
    }

    scalar dist = 0.0;
    scalar radius = 100.0;
    scalar fov = 40.0;

    if (nd->first_attribute("dist")) {
      std::string attribute(nd->first_attribute("dist")->value());
      if (!stringutils::extractFromString(attribute, dist)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of dist attribute for camera. "
                     "Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("radius")) {
      std::string attribute(nd->first_attribute("radius")->value());
      if (!stringutils::extractFromString(attribute, radius)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of radius attribute for camera. "
                     "Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("fov")) {
      std::string attribute(nd->first_attribute("fov")->value());
      if (!stringutils::extractFromString(attribute, fov)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fov attribute for camera. "
                     "Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    camera.rotation_ = rotation;
    camera.center_ = center;
    camera.dist_ = dist;
    camera.radius_ = radius;
    camera.fov_ = fov;

    return true;
  }

  return false;
}

void TwoDSceneXMLParser::loadScripts(
    rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene) {
  for (rapidxml::xml_node<>* nd = node->first_node("script"); nd;
       nd = nd->next_sibling("script")) {
    auto scr = std::make_shared<Script>();
    rapidxml::xml_attribute<>* typend = NULL;
    typend = nd->first_attribute("type");
    if (typend) {
      std::string handlertype(typend->value());
      if (handlertype == "twist")
        scr->type = Script::TWIST;
      else {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Invalid script 'type' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Invalid script 'type' attribute specified. Exiting."
                << std::endl;
      exit(1);
    }

    scr->func = Script::CUBIC;
    typend = nd->first_attribute("func");
    if (typend) {
      std::string handlertype(typend->value());
      if (handlertype == "cubic")
        scr->func = Script::CUBIC;
      else if (handlertype == "cosine")
        scr->func = Script::COSINE;
      else if (handlertype == "weno")
        scr->func = Script::WENO;
      else {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Invalid script 'func' attribute specified. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    scr->base_pos = 0.0;

    if (nd->first_attribute("x")) {
      std::string attribute(nd->first_attribute("x")->value());
      if (!stringutils::extractFromString(attribute, scr->v(0))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of x attribute for script. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("y")) {
      std::string attribute(nd->first_attribute("y")->value());
      if (!stringutils::extractFromString(attribute, scr->v(1))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of y attribute for script. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("z")) {
      std::string attribute(nd->first_attribute("z")->value());
      if (!stringutils::extractFromString(attribute, scr->v(2))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of z attribute for script. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("w")) {
      std::string attribute(nd->first_attribute("w")->value());
      if (!stringutils::extractFromString(attribute, scr->v(3))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of w attribute for script. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr->transform_with_origin = false;
    if (nd->first_attribute("ox")) {
      std::string attribute(nd->first_attribute("ox")->value());
      if (!stringutils::extractFromString(attribute, scr->origin(0))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of ox attribute for script. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
      scr->transform_with_origin = true;
    }
    if (nd->first_attribute("oy")) {
      std::string attribute(nd->first_attribute("oy")->value());
      if (!stringutils::extractFromString(attribute, scr->origin(1))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of oy attribute for script. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
      scr->transform_with_origin = true;
    }
    if (nd->first_attribute("oz")) {
      std::string attribute(nd->first_attribute("oz")->value());
      if (!stringutils::extractFromString(attribute, scr->origin(2))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of oz attribute for script. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
      scr->transform_with_origin = true;
    }

    if (nd->first_attribute("start")) {
      std::string attribute(nd->first_attribute("start")->value());
      if (!stringutils::extractFromString(attribute, scr->start)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of start attribute for script. "
                     "Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("end")) {
      std::string attribute(nd->first_attribute("end")->value());
      if (!stringutils::extractFromString(attribute, scr->end)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of end attribute for script. "
                     "Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr->ease_start = scr->ease_end = (scr->end - scr->start) / 3.0;

    if (nd->first_attribute("easestart")) {
      std::string attribute(nd->first_attribute("easestart")->value());
      if (!stringutils::extractFromString(attribute, scr->ease_start)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of start attribute for script. "
                     "Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    if (nd->first_attribute("easeend")) {
      std::string attribute(nd->first_attribute("easeend")->value());
      if (!stringutils::extractFromString(attribute, scr->ease_end)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of start attribute for script. "
                     "Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr->amplitude = 1.0;
    if (nd->first_attribute("amplitude")) {
      std::string attribute(nd->first_attribute("amplitude")->value());
      if (!stringutils::extractFromString(attribute, scr->amplitude)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of amplitude attribute for "
                     "script. Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr->base_dt = 1.0;
    if (nd->first_attribute("dt")) {
      std::string attribute(nd->first_attribute("dt")->value());
      if (!stringutils::extractFromString(attribute, scr->base_dt)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of dt attribute for script. Value "
                     "must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if (nd->first_attribute("base")) {
      std::string attribute(nd->first_attribute("base")->value());
      std::vector<std::string> bases = stringutils::split(attribute, ' ');

      scr->base_vertices.reserve(bases.size());
      for (const std::string& str : bases) {
        scalar y = 0;
        stringutils::extractFromString(str, y);
        scr->base_vertices.push_back(y);
      }
    }

    scr->frequency = 1.0;
    if (nd->first_attribute("frequency")) {
      std::string attribute(nd->first_attribute("frequency")->value());
      if (!stringutils::extractFromString(attribute, scr->frequency)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of frequency attribute for "
                     "script. Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr->group_index = 0;
    if (nd->first_attribute("group")) {
      std::string attribute(nd->first_attribute("group")->value());
      if (!stringutils::extractFromString(attribute, scr->group_index)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of i attribute for script. Value "
                     "must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr->transform_global = false;
    if (nd->first_attribute("global")) {
      std::string attribute(nd->first_attribute("global")->value());
      if (!stringutils::extractFromString(attribute, scr->transform_global)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of global attribute for script. "
                     "Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scr->script_node = 0;
    if (nd->first_attribute("node")) {
      std::string attribute(nd->first_attribute("node")->value());
      if (!stringutils::extractFromString(attribute, scr->script_node)) {
        std::cerr << outputmod::startred
          << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
          << " Failed to parse value of node attribute for script. Value "
          "must be integer. Exiting."
          << std::endl;
        exit(1);
      }
    }

    scr->strand = 0;
    if (nd->first_attribute("strand")) {
      std::string attribute(nd->first_attribute("strand")->value());
      if (!stringutils::extractFromString(attribute, scr->strand)) {
        std::cerr << outputmod::startred
          << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
          << " Failed to parse value of strand attribute for script. Value "
          "must be integer. Exiting."
          << std::endl;
        exit(1);
      }
    }

    scr->m_scene = twodscene;
    twodscene->insertScript(scr);
  }
}

void TwoDSceneXMLParser::loadXMLFile(const std::string& filename,
                                     std::vector<char>& xmlchars,
                                     rapidxml::xml_document<>& doc) {
  // Attempt to read the text from the user-specified xml file
  std::string filecontents;
  if (!loadTextFileIntoString(filename, filecontents)) {
    std::cerr << outputmod::startred
              << "ERROR IN TWODSCENEXMLPARSER:" << outputmod::endred
              << " XML scene file " << filename << ". Failed to read file."
              << std::endl;
    exit(1);
  }

  // Copy string into an array of characters for the xml parser
  for (int i = 0; i < (int)filecontents.size(); ++i)
    xmlchars.push_back(filecontents[i]);
  xmlchars.push_back('\0');

  // Initialize the xml parser with the character vector
  doc.parse<0>(&xmlchars[0]);
}

bool TwoDSceneXMLParser::loadTextFileIntoString(const std::string& filename,
                                                std::string& filecontents) {
  // Attempt to open the text file for reading
  std::ifstream textfile(filename.c_str(), std::ifstream::in);
  if (!textfile) return false;

  // Read the entire file into a single string
  std::string line;
  while (getline(textfile, line)) filecontents.append(line);

  textfile.close();

  return true;
}

void TwoDSceneXMLParser::loadSimulationType(rapidxml::xml_node<>* node,
                                            std::string& simtype) {
  assert(node != NULL);
  rapidxml::xml_node<>* nd = node->first_node("simtype");

  if (node->first_node("simtype"))
    if (node->first_node("simtype")->first_attribute("type"))
      simtype = node->first_node("simtype")->first_attribute("type")->value();
}

void TwoDSceneXMLParser::loadHairPose(
    rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene) {
  for (rapidxml::xml_node<>* nd = node->first_node("pose"); nd;
       nd = nd->next_sibling("pose")) {
    if (nd->first_attribute("path")) {
      std::string file_path(nd->first_attribute("path")->value());

      std::ifstream ifs(file_path.c_str());

      if (ifs.fail()) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to read file: " << file_path << ". Exiting."
                  << std::endl;
        exit(1);
      }

      std::string line;
      bool start_data = false;
      std::vector<Vector4s> pos;

      int idx_x = 0;
      int idx_y = 1;
      int idx_z = 2;
      int idx_theta = 3;
      int idx_seg = 4;
      int idx_ra = 5;
      int idx_ha = 7;
      int idx_actual = 10;

      int prop_order = 0;

      while (getline(ifs, line)) {
        std::vector<std::string> data;
        stringutils::split(line, ' ', data);

        if (start_data) {
          if (data.size() < 8) continue;

          Vector4s p;
          scalar r, h, area;
          int seg;
          int actual;

          stringutils::extractFromString(data[idx_actual], actual);
          if (!actual) continue;

          stringutils::extractFromString(data[idx_x], p(0));
          stringutils::extractFromString(data[idx_y], p(1));
          stringutils::extractFromString(data[idx_z], p(2));
          stringutils::extractFromString(data[idx_theta], p(3));
          stringutils::extractFromString(data[idx_seg], seg);
          stringutils::extractFromString(data[idx_ra], r);
          stringutils::extractFromString(data[idx_ha], h);

          h -= r;
          area = M_PI * (h + 2.0 * r) * h;

          pos.push_back(p);
        } else if (line.substr(0, 10) == "end_header") {
          start_data = true;
        } else {
          if (data.size() == 0) continue;

          if (data[0] == "property") {
            if (data.size() != 3) continue;

            if (data[2] == "x")
              idx_x = prop_order;
            else if (data[2] == "y")
              idx_y = prop_order;
            else if (data[2] == "z")
              idx_z = prop_order;
            else if (data[2] == "theta")
              idx_theta = prop_order;
            else if (data[2] == "segment")
              idx_seg = prop_order;
            else if (data[2] == "ra")
              idx_ra = prop_order;
            else if (data[2] == "ha")
              idx_ha = prop_order;
            else if (data[2] == "actual")
              idx_actual = prop_order;

            prop_order++;
          } else if (data[0] == "comment") {
            continue;
          }
        }
      }

      const int total_num_verts = pos.size();
      for (int i = 0; i < total_num_verts; ++i) {
        twodscene->setPosition(i, pos[i].segment<3>(0));
        twodscene->setTwist(i, pos[i](3));
      }

      ifs.close();
    } else {
      continue;
    }
  }

  twodscene->updateRestPos();
}

void TwoDSceneXMLParser::loadHairs(rapidxml::xml_node<>* node,
                                   const std::shared_ptr<TwoDScene>& twodscene,
                                   const scalar& dt) {
  assert(node != NULL);
  // Count the number of particles, edges, and strands
  int numstrands = 0;
  int numedges = twodscene->getNumEdges();
  int numparticles = twodscene->getNumParticles();

  for (rapidxml::xml_node<>* nd = node->first_node("hair"); nd;
       nd = nd->next_sibling("hair")) {
    int paramsIndex = -1;
    if (nd->first_attribute("params")) {
      std::string attribute(nd->first_attribute("params")->value());
      if (!stringutils::extractFromString(attribute, paramsIndex)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of params (ElasticParameters "
                     "index) for hair "
                  << numstrands << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of params (ElasticParameters index) "
                   "for hair "
                << numstrands << ". Exiting." << std::endl;
      exit(1);
    }

    if (paramsIndex == -1) continue;
    auto& params = twodscene->getElasticParameters(paramsIndex);

    int start = 0;
    if (nd->first_attribute("start")) {
      std::string attribute(nd->first_attribute("start")->value());
      if (!stringutils::extractFromString(attribute, start)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of start attribute for hair "
                  << numstrands << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    int count = 0;
    if (nd->first_attribute("count")) {
      std::string attribute(nd->first_attribute("count")->value());
      if (!stringutils::extractFromString(attribute, count)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of count attribute for hair "
                  << numstrands << ". Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    std::vector<int> particle_indices;
    std::vector<int> edge_indices;
    VecX particle_radius;

    if (count == 0) {
      for (rapidxml::xml_node<>* subnd = nd->first_node("p"); subnd;
           subnd = subnd->next_sibling("p")) {
        int id = -1;
        if (subnd->first_attribute("i")) {
          std::string attribute(subnd->first_attribute("i")->value());
          if (!stringutils::extractFromString(attribute, id)) {
            std::cerr << outputmod::startred
                      << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                      << " Failed to parse value of count attribute for p id "
                      << numstrands << ". Value must be integer. Exiting."
                      << std::endl;
            exit(1);
          }
        } else {
          std::cerr << outputmod::startred
                    << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                    << " Failed to parse value of count attribute for p id "
                    << numstrands << ". No attribute id." << std::endl;
          exit(1);
        }
        particle_indices.push_back(id);
      }

      count = (int)particle_indices.size();
      int num_newedges = count - 1;

      twodscene->conservativeResizeEdges(numedges + num_newedges);

      particle_radius.resize(count * 2);

      const VectorXs& scene_radius = twodscene->getRadius();

      for (int i = 0; i < count; ++i) {
        const int pidx = particle_indices[i];

        if (scene_radius(pidx * 2 + 0) == 0.0 ||
            scene_radius(pidx * 2 + 1) == 0.0) {
          const scalar radius_A = params->getRadiusA(0);
          const scalar radius_B = params->getRadiusB(0);
          twodscene->setRadius(pidx, radius_A, radius_B);
        }

        particle_radius.segment<2>(i * 2) = scene_radius.segment<2>(pidx * 2);
      }

      for (int i = 1; i < count; ++i) {
        std::pair<int, int> newedge(particle_indices[i - 1],
                                    particle_indices[i]);
        twodscene->setEdge(numedges, newedge);
        twodscene->setEdgeRestLength(numedges,
                                     (twodscene->getPosition(newedge.first) -
                                      twodscene->getPosition(newedge.second))
                                         .norm());
        // twodscene->setEdgeToParameter( numedges, paramsIndex );
        edge_indices.push_back(numedges);
        ++numedges;
      }

      for (int i = 0; i < count; ++i) {
        const int pidx = particle_indices[i];

        if (i == count - 1) {
          twodscene->setTipVerts(pidx, true);
        } else {
          twodscene->setTipVerts(pidx, false);
        }

        const scalar radius_A = particle_radius(i * 2 + 0);
        const scalar radius_B = particle_radius(i * 2 + 1);

        const scalar original_vol = twodscene->getVol()(particle_indices[i]);
        scalar vol = twodscene->getParticleRestLength(particle_indices[i]) *
                     M_PI * radius_A * radius_B;
        twodscene->setVolume(particle_indices[i], original_vol + vol);
        const scalar original_mass = twodscene->getM()(particle_indices[i] * 4);
        scalar mass = params->m_density * vol * params->m_restVolumeFraction;
        twodscene->setMass(
            particle_indices[i], original_mass + mass,
            0.25 * mass * (radius_A * radius_A + radius_B * radius_B));
        twodscene->setTwist(particle_indices[i],
                            twodscene->getSimInfo().use_twist);

        twodscene->setVolumeFraction(particle_indices[i],
          params->m_restVolumeFraction);
      }

    } else {
      int num_newedges = count - 1;

      twodscene->conservativeResizeEdges(numedges + num_newedges);

      particle_radius.resize(count * 2);

      const VectorXs& scene_radius = twodscene->getRadius();

      for (int i = 0; i < count; ++i) {
        const int pidx = start + i;

        if (scene_radius(pidx * 2 + 0) == 0.0 ||
            scene_radius(pidx * 2 + 1) == 0.0) {
          const scalar radius_A = params->getRadiusA(0);
          const scalar radius_B = params->getRadiusB(0);
          twodscene->setRadius(pidx, radius_A, radius_B);
        }

        particle_radius.segment<2>(i * 2) = scene_radius.segment<2>(pidx * 2);
      }

      for (int i = 1; i < count; ++i) {
        const int vtx = start + i;
        std::pair<int, int> newedge(vtx - 1, vtx);
        twodscene->setEdge(numedges, newedge);
        twodscene->setEdgeRestLength(numedges,
                                     (twodscene->getPosition(newedge.first) -
                                      twodscene->getPosition(newedge.second))
                                         .norm());
        // twodscene->setEdgeToParameter( numedges, paramsIndex );
        edge_indices.push_back(numedges);
        ++numedges;
      }

      particle_indices.resize(count);

      for (int i = 0; i < count; ++i) {
        const int pidx = start + i;

        if (i == count - 1) {
          twodscene->setTipVerts(pidx, true);
        } else {
          twodscene->setTipVerts(pidx, false);
        }

        particle_indices[i] = pidx;

        const scalar radius_A = particle_radius(i * 2 + 0);
        const scalar radius_B = particle_radius(i * 2 + 1);

        const scalar original_vol = twodscene->getVol()(particle_indices[i]);
        scalar vol = twodscene->getParticleRestLength(particle_indices[i]) *
                     M_PI * radius_A * radius_B;
        twodscene->setVolume(particle_indices[i], original_vol + vol);
        const scalar original_mass = twodscene->getM()(particle_indices[i] * 4);
        scalar mass = params->m_density * vol * params->m_restVolumeFraction;
        twodscene->setMass(
            particle_indices[i], original_mass + mass,
            0.25 * mass * (radius_A * radius_A + radius_B * radius_B));
        twodscene->setTwist(particle_indices[i],
                            twodscene->getSimInfo().use_twist);

        twodscene->setVolumeFraction(particle_indices[i],
          params->m_restVolumeFraction);
      }
    }

    // instance new ElasticParameters
    const int idx_sp = twodscene->getNumElasticParameters();
    twodscene->insertElasticParameters(std::make_shared<ElasticParameters>(
        particle_radius, params->m_youngsModulus.get(),
        params->m_shearModulus.get(), params->m_stretchingMultiplier,
        params->m_collisionMultiplier, params->m_attachMultiplier,
        params->m_density, params->m_viscosity, params->m_baseRotation.get(),
        dt, params->m_friction_alpha, params->m_friction_beta,
        params->m_restVolumeFraction, params->m_accumulateWithViscous,
        params->m_accumulateViscousOnlyForBendingModes,
        params->m_postProjectFixed, params->m_useApproxJacobian,
        params->m_useTournierJacobian, params->m_straightHairs,
        params->m_color));

    for (int eidx : edge_indices) {
      twodscene->setEdgeToParameter(eidx, idx_sp);
    }

    auto strandf = std::make_shared<StrandForce>(
      twodscene, particle_indices, idx_sp, numstrands);

    twodscene->insertForce(strandf);
    twodscene->insertStrandForce(strandf);

    VectorXi solve_group(particle_indices.size());
    for (int i = 0; i < (int)particle_indices.size(); ++i)
      solve_group(i) = particle_indices[i];

    twodscene->insertSolveGroup(solve_group);

    ++numstrands;
  }
}

void TwoDSceneXMLParser::loadSimInfo(
    rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene) {
  SimInfo info;
  info.yazdchi_power = 1.6;
  info.rest_volume_fraction = 0.4;
  info.lambda = 2.0;
  info.elasto_flip_asym_coeff = 0.996;
  info.elasto_flip_coeff = 0.0;
  info.elasto_advect_coeff = 0.996;
  info.solve_solid = true; 
  info.use_twist = true;
  info.use_amgpcg_solid = false;
  info.use_pcr = true;
  info.use_lagrangian_mpm = false;
  info.use_cosolve_angular = false;
  info.iteration_print_step = 0;

  twodscene->setSimInfo(info);
}

void TwoDSceneXMLParser::loadElasticParameters(
    rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene,
    const scalar& dt) {
  /* Available options, if not defined take default below
   <ElasticParameters>
   <radius value="" />
   <youngsModulus value="" />
   <shearModulus value="" />
   <density value="" />
   <viscosity value="" />
   <baseRotation value=""/>
   <accumulateWithViscous value="1" />
   <accumulateViscousOnlyForBendingModes value="1" />
   <variableRadiusHair value="1" />
   <straightHairs value="0" />
   <friction_angle value="0"/>
   </ElasticParameters>
   */

  rapidxml::xml_node<>* nd;

  int paramsCount = 0;
  for (nd = node->first_node("ElasticParameters"); nd;
       nd = nd->next_sibling("ElasticParameters")) {
    // default values:
    scalar radius = 0.018;
    scalar biradius = 0.018;
    scalar restVolumeFraction = twodscene->getSimInfo().rest_volume_fraction;
    scalar YoungsModulus = 6.687e5;
    scalar shearModulus = 2.476e5;
    scalar density = 1.3;
    scalar viscosity = 0.0;
    scalar stretchingMultiplier = 1.0;
    scalar collisionMultiplier = 1e-3;
    scalar attachMultiplier = 1e-3;
    scalar baseRotation = 0.;
    bool accumulateWithViscous = false;
    bool accumulateViscousOnlyForBendingModes = true;
    bool postProjectFixed = false;
    bool useApproxJacobian = false;
    bool useTournierJacobian = true;
    scalar straightHairs = 1.;
    Vec3 haircolor = Vec3(0, 0, 0);
    rapidxml::xml_node<>* subnd;

    if ((subnd = nd->first_node("haircolor"))) {
      std::string attributer(subnd->first_attribute("r")->value());
      if (!stringutils::extractFromString(attributer, haircolor(0))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of haircolor attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
      std::string attributeg(subnd->first_attribute("g")->value());
      if (!stringutils::extractFromString(attributeg, haircolor(1))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of haircolor attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
      std::string attributeb(subnd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attributeb, haircolor(2))) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of haircolor attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("radius"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, radius)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of radius attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("biradius"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, biradius)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of biradius attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("restVolumeFraction"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, restVolumeFraction)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of poreRadius attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("youngsModulus"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, YoungsModulus)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of youngsModulus attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("stretchingMultiplier"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, stretchingMultiplier)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of stretchingMultiplier attribute "
                     "for ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("collisionMultiplier"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, collisionMultiplier)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of collisionMultiplier attribute "
                     "for ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("attachMultiplier"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, attachMultiplier)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of attachMultiplier attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scalar poissonRatio = 0.35;
    if ((subnd = nd->first_node("poissonRatio"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, poissonRatio)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of poissonRatio attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    shearModulus = YoungsModulus / ((1.0 + poissonRatio) * 2.0);

    if ((subnd = nd->first_node("density"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, density)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of density attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("viscosity"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, viscosity)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of viscosity attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("baseRotation"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, baseRotation)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of baseRotation attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("accumulateWithViscous"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, accumulateWithViscous)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of accumulateWithViscous "
                     "attribute for ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("accumulateViscousOnlyForBendingModes"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(
              attribute, accumulateViscousOnlyForBendingModes)) {
        std::cerr
            << outputmod::startred
            << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
            << " Failed to parse value of accumulateViscousOnlyForBendingModes "
               "attribute for ElasticParameters "
            << paramsCount << ". Value must be numeric. Exiting." << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("postProjectFixed"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, postProjectFixed)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of postProjectFixed attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("useApproxJacobian"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, useApproxJacobian)) {
        std::cerr
            << outputmod::startred
            << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
            << " Failed to parse value of useApproxJacobian attribute for "
               "ElasticParameters "
            << paramsCount << ". Value must be boolean. Exiting." << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("useTournierJacobian"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, useTournierJacobian)) {
        std::cerr
            << outputmod::startred
            << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
            << " Failed to parse value of useTournierJacobian attribute for "
               "ElasticParameters "
            << paramsCount << ". Value must be boolean. Exiting." << std::endl;
        exit(1);
      }
    }

    if ((subnd = nd->first_node("straightHairs"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, straightHairs)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of straightHairs attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scalar friction_angle = 0;
    if ((subnd = nd->first_node("frictionAngle"))) {
      std::string attribute(subnd->first_attribute("value")->value());
      if (!stringutils::extractFromString(attribute, friction_angle)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of friction alpha attribute for "
                     "ElasticParameters "
                  << paramsCount << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    friction_angle =
        std::max(0., std::min(90.0, friction_angle)) / 180.0 * M_PI;
    const scalar friction_alpha =
        1.6329931619 * sin(friction_angle) / (3.0 - sin(friction_angle));
    const scalar friction_beta = tan(friction_angle);

    VecX rad_vec(2);
    rad_vec(0) = radius;
    rad_vec(1) = biradius;

    twodscene->insertElasticParameters(std::make_shared<ElasticParameters>(
        rad_vec, YoungsModulus, shearModulus, stretchingMultiplier,
        collisionMultiplier, attachMultiplier, density, viscosity, baseRotation,
        dt, friction_alpha, friction_beta, restVolumeFraction,
        accumulateWithViscous, accumulateViscousOnlyForBendingModes,
        postProjectFixed, useApproxJacobian, useTournierJacobian, straightHairs,
        haircolor));
    ++paramsCount;
  }
}

void TwoDSceneXMLParser::loadParticles(
    rapidxml::xml_node<>* node, const std::shared_ptr<TwoDScene>& twodscene,
    int& maxgroup) {
  // Count the number of particles
  int old_num_particles = twodscene->getNumParticles();
  int numparticles = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("particle"); nd;
       nd = nd->next_sibling("particle"))
    ++numparticles;

  twodscene->resizeParticleSystem(numparticles);

  // std::cout << "Num particles " << numparticles << std::endl;

  maxgroup = 0;

  int particle = 0;
  for (rapidxml::xml_node<>* nd = node->first_node("particle"); nd;
       nd = nd->next_sibling("particle")) {
    // Extract the particle's initial position
    Vector3s pos = Vector3s::Zero();
    if (nd->first_attribute("x")) {
      std::string position(nd->first_attribute("x")->value());
      if (!stringutils::readList(position, ' ', pos)) {
        std::cerr << "Failed to load x, y, and z positions for particle "
                  << particle << std::endl;
        exit(1);
      }
    } else {
      std::cerr
          << "Failed to find x, y, and z position attributes for particle "
          << particle << std::endl;
      exit(1);
    }
    twodscene->setPosition(particle, pos);

    // Extract the particle's initial velocity
    Vector3s vel = Vector3s::Zero();
    if (nd->first_attribute("v")) {
      std::string velocity(nd->first_attribute("v")->value());
      if (!stringutils::readList(velocity, ' ', vel)) {
        std::cerr << "Failed to load x, y, and z velocities for particle "
                  << particle << std::endl;
        exit(1);
      }
    }

    twodscene->setVelocity(particle, vel);

    // parse theta
    scalar theta = 0.0;
    if (nd->first_attribute("theta")) {
      std::string attribute(nd->first_attribute("theta")->value());
      if (!stringutils::extractFromString(attribute, theta)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fixed attribute for particle "
                  << particle << ". Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    twodscene->setTheta(particle, theta);

    scalar omega = 0.0;
    if (nd->first_attribute("omega")) {
      std::string attribute(nd->first_attribute("omega")->value());
      if (!stringutils::extractFromString(attribute, omega)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fixed attribute for particle "
                  << particle << ". Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    twodscene->setOmega(particle, omega);

    // Determine if the particle is fixed
    int fixed = 0;
    if (nd->first_attribute("fixed")) {
      std::string attribute(nd->first_attribute("fixed")->value());
      if (!stringutils::extractFromString(attribute, fixed)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fixed attribute for particle "
                  << particle << ". Value must be boolean. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    twodscene->setFixed(particle, (unsigned char)(fixed & 0xFFU));
    twodscene->setTwist(particle, false);

    // Extract the particle's radius, if present
    scalar radius = 0.0;
    if (nd->first_attribute("radius")) {
      std::string attribute(nd->first_attribute("radius")->value());
      if (!stringutils::extractFromString(attribute, radius)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse radius attribute for particle "
                  << particle << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scalar biradius = radius;
    if (nd->first_attribute("biradius")) {
      std::string attribute(nd->first_attribute("biradius")->value());
      if (!stringutils::extractFromString(attribute, biradius)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse biradius attribute for particle "
                  << particle << ". Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    twodscene->setRadius(particle, radius, biradius);

    scalar vol = 0.0;
    if (nd->first_attribute("vol")) {
      std::string attribute(nd->first_attribute("vol")->value());
      if (!stringutils::extractFromString(attribute, vol)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse vol attribute for particle " << particle
                  << ". Value must be scalar. Exiting." << std::endl;
        exit(1);
      }
    }
    twodscene->setVolume(particle, vol);

    scalar fvol = 0.0;
    if (nd->first_attribute("fvol")) {
      std::string attribute(nd->first_attribute("fvol")->value());
      if (!stringutils::extractFromString(attribute, fvol)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse fvol attribute for particle " << particle
                  << ". Value must be scalar. Exiting." << std::endl;
        exit(1);
      }
    }
    twodscene->setFluidVolume(particle, fvol);

    // Extract the particle's mass
    scalar mass = 0.0;
    if (nd->first_attribute("m")) {
      std::string attribute(nd->first_attribute("m")->value());
      if (!stringutils::extractFromString(attribute, mass)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of m attribute for particle "
                  << particle << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    twodscene->setMass(particle, mass, 0.0);

    scalar fmass = 0.0;
    if (nd->first_attribute("fm")) {
      std::string attribute(nd->first_attribute("fm")->value());
      if (!stringutils::extractFromString(attribute, fmass)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of fm attribute for particle "
                  << particle << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    twodscene->setFluidMass(particle, fmass, 0.0);

    if (nd->first_attribute("state")) {
      std::string attribute(nd->first_attribute("state")->value());
      if (attribute == "liquid") {
        twodscene->getFluidIndices().push_back(particle);
      }
    }

    scalar vf = 1.0;
    if (nd->first_attribute("vf")) {
      std::string attribute(nd->first_attribute("vf")->value());
      if (!stringutils::extractFromString(attribute, vf)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of vf attribute for particle "
                  << particle << ". Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }
    twodscene->setVolumeFraction(particle, vf);

    // std::cout << "Particle: " << particle << "    x: " << pos.transpose() <<
    // "   v: " << vel.transpose() << "   m: " << mass << "   fixed: " << fixed
    // << std::endl; std::cout << tags[particle] << std::endl;
    ++particle;
  }
}

void TwoDSceneXMLParser::loadSceneTag(rapidxml::xml_node<>* node,
                                      std::string& scenetag) {
  assert(node != NULL);

  if (node->first_node("scenetag")) {
    if (node->first_node("scenetag")->first_attribute("tag")) {
      scenetag = node->first_node("scenetag")->first_attribute("tag")->value();
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of tag attribute for scenetag. "
                   "Value must be string. Exiting."
                << std::endl;
      exit(1);
    }
  }
}

void TwoDSceneXMLParser::loadIntegrator(
    rapidxml::xml_node<>* node, std::shared_ptr<SceneStepper>& scenestepper,
    scalar& dt) {
  assert(node != NULL);

  dt = -1.0;

  // Attempt to locate the integrator node
  rapidxml::xml_node<>* nd = node->first_node("integrator");
  if (nd == NULL) {
    std::cerr << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " No integrator specified. Exiting." << std::endl;
    exit(1);
  }

  // Attempt to load the integrator type
  rapidxml::xml_attribute<>* typend = nd->first_attribute("type");
  if (typend == NULL) {
    std::cerr << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " No integrator 'type' attribute specified. Exiting."
              << std::endl;
    exit(1);
  }
  std::string integratortype(typend->value());

  if (integratortype == "linearized-implicit-euler") {
    rapidxml::xml_attribute<>* subnd;

    int manifoldsubsteps = 4;
    subnd = nd->first_attribute("manifoldsubsteps");
    if (subnd) {
      if (!stringutils::extractFromString(std::string(subnd->value()),
                                          manifoldsubsteps)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse 'manifoldsubsteps' attribute for "
                     "integrator. Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scalar criterion = 1e-6;
    subnd = nd->first_attribute("criterion");
    if (subnd) {
      if (!stringutils::extractFromString(std::string(subnd->value()),
                                          criterion)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse 'criterion' attribute for integrator. "
                     "Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scalar pressure_criterion = criterion;
    subnd = nd->first_attribute("pressurecriterion");
    if (subnd) {
      if (!stringutils::extractFromString(std::string(subnd->value()),
                                          pressure_criterion)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse 'pressurecriterion' attribute for "
                     "integrator. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scalar viscous_criterion = criterion;
    subnd = nd->first_attribute("viscouscriterion");
    if (subnd) {
      if (!stringutils::extractFromString(std::string(subnd->value()),
                                          viscous_criterion)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse 'pressurecriterion' attribute for "
                     "integrator. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scalar quasi_static_criterion = criterion;
    subnd = nd->first_attribute("quasistaticcriterion");
    if (subnd) {
      if (!stringutils::extractFromString(std::string(subnd->value()),
                                          quasi_static_criterion)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse 'quasistaticcriterion' attribute for "
                     "integrator. Value must be numeric. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    int maxiters = 100;
    subnd = nd->first_attribute("maxiters");
    if (subnd) {
      if (!stringutils::extractFromString(std::string(subnd->value()),
                                          maxiters)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse 'maxiters' attribute for integrator. "
                     "Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    int viscositysubsteps = 1;
    subnd = nd->first_attribute("viscositysubsteps");
    if (subnd) {
      if (!stringutils::extractFromString(std::string(subnd->value()),
                                          viscositysubsteps)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse 'viscositysubsteps' attribute for "
                     "integrator. Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    int surftensionsubsteps = 1;
    subnd = nd->first_attribute("surftensionsubsteps");
    if (subnd) {
      if (!stringutils::extractFromString(std::string(subnd->value()),
                                          surftensionsubsteps)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse 'surftensionsubsteps' attribute for "
                     "integrator. Value must be integer. Exiting."
                  << std::endl;
        exit(1);
      }
    }

    scenestepper = std::make_shared<LinearizedImplicitEuler>(
        criterion, pressure_criterion, quasi_static_criterion,
        viscous_criterion, maxiters, manifoldsubsteps, viscositysubsteps,
        surftensionsubsteps);
  } else {
    std::cerr << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Invalid integrator 'type' attribute specified. Exiting."
              << std::endl;
    exit(1);
  }

  // Attempt to load the timestep
  rapidxml::xml_attribute<>* dtnd = nd->first_attribute("dt");
  if (dtnd == NULL) {
    std::cerr << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " No integrator 'dt' attribute specified. Exiting."
              << std::endl;
    exit(1);
  }

  dt = std::numeric_limits<scalar>::signaling_NaN();
  if (!stringutils::extractFromString(std::string(dtnd->value()), dt)) {
    std::cerr << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse 'dt' attribute for integrator. Value must "
                 "be numeric. Exiting."
              << std::endl;
    exit(1);
  }

  bool useApic = false;
  rapidxml::xml_attribute<>* apnd = nd->first_attribute("apic");
  if (apnd) {
    if (!stringutils::extractFromString(std::string(apnd->value()), useApic)) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse 'apic' attribute for integrator. Value "
                   "must be bool. Exiting."
                << std::endl;
      exit(1);
    }
  }
  scenestepper->setUseApic(useApic);

  // std::cout << "Integrator: " << (*scenestepper)->getName() << "   dt: " <<
  // dt << std::endl;
}

void TwoDSceneXMLParser::loadMaxTime(rapidxml::xml_node<>* node,
                                     scalar& max_t) {
  assert(node != NULL);

  // Attempt to locate the duraiton node
  rapidxml::xml_node<>* nd = node->first_node("duration");
  if (nd == NULL) {
    std::cerr << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " No duration specified. Exiting." << std::endl;
    exit(1);
  }

  // Attempt to load the duration value
  rapidxml::xml_attribute<>* timend = nd->first_attribute("time");
  if (timend == NULL) {
    std::cerr << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " No duration 'time' attribute specified. Exiting."
              << std::endl;
    exit(1);
  }

  max_t = std::numeric_limits<scalar>::signaling_NaN();
  if (!stringutils::extractFromString(std::string(timend->value()), max_t)) {
    std::cerr << outputmod::startred
              << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
              << " Failed to parse 'time' attribute for duration. Value must "
                 "be numeric. Exiting."
              << std::endl;
    exit(1);
  }
}

void TwoDSceneXMLParser::loadViewport(rapidxml::xml_node<>* node,
                                      renderingutils::Viewport& view) {
  assert(node != NULL);

  if (node->first_node("viewport")) {
    rapidxml::xml_attribute<>* cx =
        node->first_node("viewport")->first_attribute("cx");
    if (cx == NULL) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " No viewport 'cx' attribute specified. Exiting."
                << std::endl;
      exit(1);
    }
    if (!stringutils::extractFromString(std::string(cx->value()), view.cx)) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse 'cx' attribute for viewport. Value must "
                   "be scalar. Exiting."
                << std::endl;
      exit(1);
    }
    rapidxml::xml_attribute<>* cy =
        node->first_node("viewport")->first_attribute("cy");
    if (cy == NULL) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " No viewport 'cy' attribute specified. Exiting."
                << std::endl;
      exit(1);
    }
    if (!stringutils::extractFromString(std::string(cy->value()), view.cy)) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse 'cy' attribute for viewport. Value must "
                   "be scalar. Exiting."
                << std::endl;
      exit(1);
    }
    rapidxml::xml_attribute<>* size =
        node->first_node("viewport")->first_attribute("size");
    if (size == NULL) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " No viewport 'size' attribute specified. Exiting."
                << std::endl;
      exit(1);
    }
    if (!stringutils::extractFromString(std::string(size->value()),
                                        view.size)) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse 'size' attribute for viewport. Value must "
                   "be scalar. Exiting."
                << std::endl;
      exit(1);
    }
  }
}

void TwoDSceneXMLParser::loadMaxSimFrequency(rapidxml::xml_node<>* node,
                                             scalar& max_freq) {
  assert(node != NULL);

  // Attempt to locate the duraiton node
  if (node->first_node("maxsimfreq")) {
    // Attempt to load the duration value
    rapidxml::xml_attribute<>* atrbnde =
        node->first_node("maxsimfreq")->first_attribute("max");
    if (atrbnde == NULL) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " No maxsimfreq 'max' attribute specified. Exiting."
                << std::endl;
      exit(1);
    }

    if (!stringutils::extractFromString(std::string(atrbnde->value()),
                                        max_freq)) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse 'max' attribute for maxsimfreq. Value "
                   "must be scalar. Exiting."
                << std::endl;
      exit(1);
    }
  }
}

void TwoDSceneXMLParser::loadBackgroundColor(rapidxml::xml_node<>* node,
                                             renderingutils::Color& color) {
  if (rapidxml::xml_node<>* nd = node->first_node("backgroundcolor")) {
    // Read in the red color channel
    double red = -1.0;
    if (nd->first_attribute("r")) {
      std::string attribute(nd->first_attribute("r")->value());
      if (!stringutils::extractFromString(attribute, red)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of r attribute for "
                     "backgroundcolor. Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of r attribute for backgroundcolor. "
                   "Exiting."
                << std::endl;
      exit(1);
    }

    if (red < 0.0 || red > 1.0) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of r attribute for backgroundcolor. "
                   "Invalid color specified. Valid range is "
                << 0.0 << "..." << 1.0 << std::endl;
      exit(1);
    }

    // Read in the green color channel
    double green = -1.0;
    if (nd->first_attribute("g")) {
      std::string attribute(nd->first_attribute("g")->value());
      if (!stringutils::extractFromString(attribute, green)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of g attribute for "
                     "backgroundcolor. Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of g attribute for backgroundcolor. "
                   "Exiting."
                << std::endl;
      exit(1);
    }

    if (green < 0.0 || green > 1.0) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of g attribute for backgroundcolor. "
                   "Invalid color specified. Valid range is "
                << 0.0 << "..." << 1.0 << std::endl;
      exit(1);
    }

    // Read in the blue color channel
    double blue = -1.0;
    if (nd->first_attribute("b")) {
      std::string attribute(nd->first_attribute("b")->value());
      if (!stringutils::extractFromString(attribute, blue)) {
        std::cerr << outputmod::startred
                  << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                  << " Failed to parse value of b attribute for "
                     "backgroundcolor. Value must be scalar. Exiting."
                  << std::endl;
        exit(1);
      }
    } else {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of b attribute for backgroundcolor. "
                   "Exiting."
                << std::endl;
      exit(1);
    }

    if (blue < 0.0 || blue > 1.0) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " Failed to parse value of b attribute for backgroundcolor. "
                   "Invalid color specified. Valid range is "
                << 0.0 << "..." << 1.0 << std::endl;
      exit(1);
    }

    // std::cout << red << "   " << green << "   " << blue << std::endl;

    color.r = red;
    color.g = green;
    color.b = blue;
  }
}

void TwoDSceneXMLParser::loadSceneDescriptionString(
    rapidxml::xml_node<>* node, std::string& description_string) {
  assert(node != NULL);

  description_string = "No description specified.";

  // Attempt to locate the integrator node
  rapidxml::xml_node<>* nd = node->first_node("description");
  if (nd != NULL) {
    rapidxml::xml_attribute<>* typend = nd->first_attribute("text");
    if (typend == NULL) {
      std::cerr << outputmod::startred
                << "ERROR IN XMLSCENEPARSER:" << outputmod::endred
                << " No text attribute specified for description. Exiting."
                << std::endl;
      exit(1);
    }
    description_string = typend->value();
  }
}
