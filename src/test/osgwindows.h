#pragma once


#include <osgDB/ReadFile>
#include <osgUtil/Optimizer>
#include <osgViewer/Viewer>


#include <iostream>

#define TEST_DATA_PATH "../test_data"


int osgwindows(int argc, char** argv)
{
    // use an ArgumentParser object to manage the program arguments.
    osg::ArgumentParser arguments(&argc, argv);

    // read the scene from the list of file specified commandline args.
    osg::ref_ptr<osg::Node> loadedModel = osgDB::readRefNodeFiles(arguments);

    // if not loaded assume no arguments passed in, try use default mode instead.

        // Load a mesh in OFF format
    std::string fileName = TEST_DATA_PATH"/OpenSceneGraph-Data/cow.osg";

    std::cout << "fileName: " << fileName << std::endl;

    if (!loadedModel) loadedModel = osgDB::readRefNodeFile(fileName);

    // if no model has been successfully loaded report failure.
    if (!loadedModel)
    {
        std::cout << arguments.getApplicationName() << ": No data loaded" << std::endl;
        return 1;
    }

    // construct the viewer.
    osgViewer::Viewer viewer;

    int xoffset = 40;
    int yoffset = 40;

    // left window + left slave camera
    {
        osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;

        traits->x = xoffset + 0;
        traits->y = yoffset + 0;
        traits->width = 600;
        traits->height = 480;
        traits->windowDecoration = true;
        traits->doubleBuffer = true;
        traits->sharedContext = 0;
        traits->readDISPLAY();
        traits->setUndefinedScreenDetailsToDefaultScreen();

        osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());

        osg::ref_ptr<osg::Camera> camera = new osg::Camera;
        camera->setGraphicsContext(gc.get());
        camera->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));
        GLenum buffer = traits->doubleBuffer ? GL_BACK : GL_FRONT;
        camera->setDrawBuffer(buffer);
        camera->setReadBuffer(buffer);

        // add this slave camera to the viewer, with a shift left of the projection matrix
        viewer.addSlave(camera.get(), osg::Matrixd::translate(1.0, 0.0, 0.0), osg::Matrixd());
    }

    // right window + right slave camera
    {
        osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
        traits->x = xoffset + 600;
        traits->y = yoffset + 0;
        traits->width = 600;
        traits->height = 480;
        traits->windowDecoration = true;
        traits->doubleBuffer = true;
        traits->sharedContext = 0;
        traits->readDISPLAY();
        traits->setUndefinedScreenDetailsToDefaultScreen();

        osg::ref_ptr<osg::GraphicsContext> gc = osg::GraphicsContext::createGraphicsContext(traits.get());

        osg::ref_ptr<osg::Camera> camera = new osg::Camera;
        camera->setGraphicsContext(gc.get());
        camera->setViewport(new osg::Viewport(0, 0, traits->width, traits->height));
        GLenum buffer = traits->doubleBuffer ? GL_BACK : GL_FRONT;
        camera->setDrawBuffer(buffer);
        camera->setReadBuffer(buffer);

        // add this slave camera to the viewer, with a shift right of the projection matrix
        viewer.addSlave(camera.get(), osg::Matrixd::translate(-1.0, 0.0, 0.0), osg::Matrixd());
    }


    // optimize the scene graph, remove redundant nodes and state etc.
    osgUtil::Optimizer optimizer;
    optimizer.optimize(loadedModel);

    // set the scene to render
    viewer.setSceneData(loadedModel);

    return viewer.run();
}