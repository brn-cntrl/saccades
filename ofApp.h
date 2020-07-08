#pragma once

#include "ofMain.h"
#include "ofxOsc.h"
#include "ofxGui.h"
#include "ofxSyphon.h"

#define PORT 31337

class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();
        void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
        void drawMesh();
        void drawNoise(float noiseF);
        void clearFbo(ofFbo buffer);
        void drawBloom(float intens);
        void drawBlurX(float blr);
        void drawBlurY(float blr);
        void drawBackbone();
        void cycleSideChains(int indexNum);
        void createNormals(ofMesh m);
        void drawToFbo(ofFbo fbo, float a, int b = 0);
        void orbitCam();
        void calculateIndices(ofMesh m, float distanceThreshold);
        //vector<string> packAminoAcids(vector<string> s);
        ofVec3f getCenter(const vector< vector<float> > &v);
        vector<string> split(string str, char delimiter);
        vector<float> getFloatVector(vector<string> &sv);
        string removeSpaces(string str);
        vector<string> getAminoAcid(ofBuffer b);
    
    //Create struct for position / acid pairing
        struct AcidPosition
        {
            vector<float> position;
            string acid;
        };
    
    //Syphon setup
    //ofxSyphonServer server;
    
    //Map for hydrophobicity values
    std::map<std::string, string> hydrovals;
    
    AcidPosition pairAcids(vector<vector<float>> pos, vector<string> acids, int i);
    
    //ofxPanel gui;
    float blurAmount = 1.8;
    float intensity = .50;
    float xVal = 0;
    float zVal = 0;
    float noiseF = 4;
    
    ofEasyCam cam;
    ofMesh mesh;
    ofMesh backbone;
    ofLight light;
    ofxOscSender sender;
    //ofxOscReceiver receiver;
    ofShader shaderBlurX;
    ofShader shaderBlurY;
    ofShader shaderNoise;
    ofShader shaderBloom;
    ofFbo fbo;
    ofFbo noiseFbo;
    ofFbo bloomFbo;
    ofFbo blurXFbo;
    ofFbo blurYFbo;
    ofNode posNode;
    ofVec3f currentLoc;
    //OSC values to send TODO: Make function and scope variables (return osc message)
    float oscX, oscY, oscZ;
    //Mesh Rotation
    float rotx, roty, rotz;
    //Distance between coordinates
    float connectionDistance = 3.7;
    //Counter to track vertices for target square
    unsigned int vertexCounter = 0;
    //Number of vertices in current PDB mesh
    int numVertices;
    //Full mesh coordinates
    vector< vector<float> > coordinates;
    //Backbone mesh coordinates
    vector< vector<float> > backboneCoords;
    //Amino acid group coordinates
    vector< vector<float> > sideChainCoords;
    //Amino acid group temp mesh
    ofMesh aagTempMesh;
    //Amino acid group meshes vector
    vector<ofMesh> sideChainGroups;
    //Counter for mesh vector
    int mvCounter = 0;
    //List of amino acids
    vector< string > aminoAcidLabels;
    //List of labels for atoms in sidechain
    vector< vector<string> > sideChainAtoms;
    //String to display amino acid
    string displayAminoAcid;
    //Temp float vector for sidechain atom coordinates
    //display Vector
    vector<string> displayVector;
    vector<float> aagFloatVector;
    //Model Rotation Speed
    float rotationSpeed = .00005f;
    //Number of side chain meshes in array
    int numSideChains;
    //Side chain ID variables
    int currentID = 0;
    int lastID = 0;
    //Hydrophobicity of current amino acid
    string currHydro;
    //tempo for sound and animation
    int tempo = 8;
    
    
    
    
};
