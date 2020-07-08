#include "ofMain.h"
#include "ofApp.h"
#include <iostream>
#include <fstream>

using namespace std;
//========================================================================
int main( ){
	ofSetupOpenGL(1920,1080,OF_WINDOW);// <-------- setup the GL context
    ofGLFWWindowSettings settings;
    settings.setGLVersion(3, 2);
	// this kicks off the running of my app
	// can be OF_WINDOW or OF_FULLSCREEN
	// pass in width and height too:
	ofRunApp(new ofApp());

}
