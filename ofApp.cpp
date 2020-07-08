/*
 Saccade / GLP1R for iMappening
 Brian Cantrell, 2019
 
 TODO: Capture at 1920 x 1080 (set display and cature in screenflow.
 Apply TV snow shader.
 
 
 */
#include "ofApp.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
using namespace std;

//------------------BEGIN SETUP---------------------------------
void ofApp::setup()
{
    //Syphon server
    //server.setName("simple server");
    
    //Setup Gui-------------------------------------------------
    /*
    gui.setup("TITLE");
    gui.add(blurAmount.set("blur",.8, .01, 3));
    gui.add(intensity.set("intensity", .5, .01, 1));
    gui.add(xVal.set("X value", 0, -600, 600));
    gui.add(zVal.set("Z value", 0, -600, 600));
     */
    
    // BEGIN: GL BOILERPLATE------------------------------------
    
    //Sync redraw with vertical refresh rate of screen
    ofSetVerticalSync(true);
    ofSetFrameRate(60);
    //Set distance of camera from object
    cam.setDistance(50);
    
    //Allocate FBOs
    fbo.allocate(1920, 1080, GL_RGBA);
    noiseFbo.allocate(1920, 1080, GL_RGBA);
    bloomFbo.allocate(1920, 1080, GL_RGBA);
    blurXFbo.allocate(2420, 1580, GL_RGBA);
    blurYFbo.allocate(2420, 1580, GL_RGBA);
    
    //clear FBO (FBO is always cleared in setup)
    
    clearFbo(fbo);
    clearFbo(noiseFbo);
    clearFbo(bloomFbo);
    clearFbo(blurXFbo);
    clearFbo(blurYFbo);
    
    //Render lighting
    light.setup();
    light.setPosition(0, 100, 250);
    //light.enable();
    
    //Set type for primary meshes
    mesh.setMode(OF_PRIMITIVE_POINTS);
    backbone.setMode(OF_PRIMITIVE_POINTS);
    //reticule.setMode(OF_PRIMITIVE_POINTS);
    
    
    //enable indexing
    mesh.enableIndices();
    backbone.enableIndices();
    //Set point size
    glPointSize(18);
    
    //Set smoothing
    glEnable(GL_POINT_SMOOTH);
    
    // END: GL BOILERPLATE--------------------------------------
    
//Preprocessor directive to load Shader for apprpriate version of OpenGL
#ifdef TARGET_OPENGLES
    shaderBlurX.load("shadersES2/shaderBlurX");
    shaderBlurY.load("shadersES2/shaderBlurY");
#else
    if(ofIsGLProgrammableRenderer())
    {
        shaderBlurX.load("shadersGL3/shaderBlurX");
        shaderBlurY.load("shadersGL3/shaderBlurY");
    }
    else
    {
        shaderBlurX.load("shadersGL2/shaderBlurX");
        shaderBlurY.load("shadersGL2/shaderBlurY");
    }
#endif
    shaderNoise.load("noise/random");
    shaderBloom.load("bloom/bloom");
    
    
    // BEGIN: STRIP OUT ATOM COORDINATES AND CONVERT TO FLOAT---
    ofBuffer file = ofBufferFromFile("positions.txt");
    ofBuffer::Lines lines = file.getLines();
    ofBuffer::Line iter = lines.begin();
    
    //Get the first ID number
    string firstIDstring = removeSpaces((*iter).substr(27, 4));
    
    //Get the amino acid type for first group
    string currentAA = (*iter).substr(20,3);

    //Push amino acid name to string vector
    aminoAcidLabels.push_back(currentAA);
    
    //Temp string and string vector to hold amino acid sidechain labels
    vector< string > tempSCAVector;
    string tempSideChainAtom;
    
    //Get the ID number for the first group
    currentID = stoi(firstIDstring);
    //Set to lastID for checking at beginning of loop
    lastID = currentID;
    
    //Boolean to check if we are in the same group on loop iteration
    bool sameGroup = true;
    //Read through each line of buffer
    while(iter != lines.end())
    {
        //Create vector for coordinates
        //search through each line for '?' and create substring
        size_t first;
        size_t second;
        string sub1;
        
        //check bounds of string before assigning variables
        if((*iter).find("?") != string::npos)
        {
            size_t first = (*iter).find("?");
            size_t second = (*iter).find("?");
            sub1 = (*iter).substr(first+2, second);
        }
        //resize to lop off final characters
        sub1.resize(24);
        
        //split each string on spaces
        vector<string> splitStrings = split(removeSpaces(sub1), ' ');
        
        //push float vector into nested vector
        coordinates.push_back(getFloatVector(splitStrings));
        
        //Create and fill backbone coordinates vector
        //substring for C-alpha positions
        
        if ((*iter).substr(14,2) == "CA" || (*iter).substr(14,2) == "C " || (*iter).substr(14,2) == "N " || (*iter).substr(14,2) == "O ")
        {
            string bbString = (*iter).substr(34,23);
            vector<string> splitStrings2 = split(removeSpaces(bbString), ' ');
            vector<float> bbFloatVector = getFloatVector(splitStrings2);
            backboneCoords.push_back(bbFloatVector);
        }
        
        else
        {
            
            //Get amino acid group number from file (TODO: Get actual name for this number.
            //Parse number as int and set number to sentry variable
            //NOTE: "25, 2" gives ATOM HETATM Water etc. "27, 4" yields AA group IDs.
            string idString = removeSpaces((*iter).substr(27, 4));
            //Cast id number as int.
            currentID = stoi(idString);
            
            if(currentID == lastID)
            {
                sameGroup = true;
            } else {
                sameGroup = false;
            }
            
            if (sameGroup)
            {
                //Get the coordinates from current line.
                string aaCoordString = (*iter).substr(34,23);
                //Split and remove spaces, push to 3 element string vector
                vector<string> aaCoordStrings = split(removeSpaces(aaCoordString), ' ');
                //Get floats and push to 3 element float vector
                aagFloatVector = getFloatVector(aaCoordStrings);
                //Push the float vectors to 3D vector
                //sideChainCoords.push_back(aagFloatVector);
                //Index ofMesh vector with counter and push xyz coordinates to Mesh
                aagTempMesh.addVertex(ofPoint(aagFloatVector[0], aagFloatVector[1], aagFloatVector[2]));
                //Get sidechain atom and push to a temporary variable
                string tempString = (*iter).substr(14,4);
                tempSideChainAtom = removeSpaces(tempString);
                //cout << (*iter).substr(14,3) << endl;
                //push atom name to vector
                tempSCAVector.push_back(tempSideChainAtom);
                
            } else {
                //Push temp sidechain atom labels vector to 2d vector
                sideChainAtoms.push_back(tempSCAVector);
                
                //Push temp mesh to mesh vector
                sideChainGroups.push_back(aagTempMesh);
                
                //Clear temp mesh
                aagTempMesh.clear();
                
                //Clear temp sidechain atom vector
                tempSCAVector.clear();
                
                //Get the current Amino acid string
                currentAA = (*iter).substr(20,3);
                
                //Push to labels vector
                aminoAcidLabels.push_back(currentAA);
                
                //Get the coordinates from current line.
                string aaCoordString = (*iter).substr(34,23);
                //Split and remove spaces, push to 3 element string vector
                vector<string> aaCoordStrings = split(removeSpaces(aaCoordString), ' ');
                //Get floats and push to 3 element float vector
                aagFloatVector = getFloatVector(aaCoordStrings);
                
                //Index ofMesh vector with counter and push xyz coordinates to Mesh
                aagTempMesh.addVertex(ofPoint(aagFloatVector[0], aagFloatVector[1], aagFloatVector[2]));
                string tempString = (*iter).substr(14,4);
                tempSideChainAtom = removeSpaces(tempString);
                //cout << (*iter).substr(14,3) << endl;
                //push atom name to vector
                tempSCAVector.push_back(tempSideChainAtom);
                lastID = currentID;
            }
            
        }
        //cout << sub2 << endl;
        
        ++iter;
    }
    //Push final mesh to vector
    sideChainGroups.push_back(aagTempMesh);
    sideChainAtoms.push_back(tempSCAVector);
    tempSCAVector.clear();
    aagTempMesh.clear();
    
    //-----END: STRIP OUT ATOM COORDINATES AND CONVERT TO FLOAT---
 
    
    //Fill list of Amino Acids
    //aminoAcids = getAminoAcid(file);

    //Fill primary mesh---------------------
    for(int i = 0; i < coordinates.size(); i++)
    {
        mesh.addVertex(ofPoint(coordinates[i][0], coordinates[i][1], coordinates[i][2]));
    }
    //Fill backbone mesh
    for(int i = 0; i < backboneCoords.size(); i++)
    {
        backbone.addVertex(ofPoint(backboneCoords[i][0], backboneCoords[i][1], backboneCoords[i][2]));
    }
    
    //-----------Calculate Indices for Drawing-------------
    
    calculateIndices(mesh, connectionDistance);
    calculateIndices(backbone, connectionDistance);
    
    numSideChains = sideChainGroups.size();
    
    for(int i = 0; i < numSideChains; i++)
    {
        calculateIndices(sideChainGroups[i], connectionDistance);
        sideChainGroups[i].setMode(OF_PRIMITIVE_POINTS);
    }
    //Calculate Normals--------------------------
    //createNormals(mesh);
    
    //Get global number of vertices in current mesh
    //TODO: Maybe strip this out?
    numVertices = mesh.getNumVertices();
    
    //Setup OSC
    sender.setup("localhost", PORT);
   
    //receiver.setup(PORT);
    //AcidPosition currentAA = pairAcids(coordinates, aminoAcids, 7);
    
    //Fill hydrovals map with amino acid names and hydrophobicity values
    //Retrieved from https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html#hydro
    
    hydrovals["PHE"] = "100";
    hydrovals["ILE"] = " 99";
    hydrovals["TRP"] = " 97";
    hydrovals["LEU"] = " 97";
    hydrovals["VAL"] = " 76";
    hydrovals["MET"] = " 74";
    hydrovals["TYR"] = " 63";
    hydrovals["CYS"] = " 49";
    hydrovals["ALA"] = " 41";
    hydrovals["THR"] = " 13";
    hydrovals["HCS"] = "  0";
    hydrovals["HIS"] = "  8";
    hydrovals["SER"] = " -5";
    hydrovals["GLN"] = "-10";
    hydrovals["ARG"] = "-14";
    hydrovals["LYS"] = "-23";
    hydrovals["ASN"] = "-28";
    hydrovals["GLU"] = "-31";
    hydrovals["PRO"] = "-46";
    hydrovals["ASP"] = "-55";
    
    
    //int hydro = hydrovals.at("GLU");
    //cout << hydro << endl;
    /*
    for(std::vector< vector< string > >::const_iterator itr = sideChainAtoms.begin(); itr != sideChainAtoms.end(); ++itr)
    {
        for(std::vector<string>::const_iterator jitr = (*itr).begin(); jitr != (*itr).end(); ++jitr)
        {
            cout << *jitr << endl;
        }
    }
    */
    
}

//-----------------------END SETUP------------------------------


//---------------------BEGIN USR FUNCTIONS----------------------

//--------------------------------------------------------------
vector<string> ofApp::getAminoAcid(ofBuffer b)
{
    vector<string> v;
    ofBuffer::Lines lines = b.getLines();
    ofBuffer::Line iter = lines.begin();
    string aa;
    while(iter != lines.end())
    {
        aa = (*iter).substr(75,3);
        v.push_back(aa);
        ++iter;
    }
    return v;
}

//--------------------------------------------------------------
ofApp::AcidPosition ofApp::pairAcids(vector<vector<float>> pos, vector<string> acids, int i)
{
    AcidPosition Pair;
    
    if (pos.size() == acids.size())
    {
        Pair.position = pos[i];
        Pair.acid = acids[i];
        return Pair;
    }
    else
    {
        cout << "Array quantities mismatch." << endl;
    }
}

//--------------------------------------------------------------
void ofApp::calculateIndices(ofMesh m, float distanceThreshold)
{
    int numVerts = m.getNumVertices();
    
    for (int a=0; a<numVerts; ++a)
    {
        ofVec3f verta = m.getVertex(a);
        for (int b=a+1; b<numVerts; ++b)
        {
            ofVec3f vertb = m.getVertex(b);
            float distance = verta.distance(vertb);
            if (distance <= distanceThreshold)
            {
                m.addIndex(a);
                m.addIndex(b);
            }
        }
    }
}

//--------------------------------------------------------------
ofVec3f ofApp::getCenter(const vector< vector<float> > &v)
{
    //Set variables to hold smallest and largest vals to
    
    float smallestX = v[0][0];
    float largestX = v[0][0];
    
    float smallestY = v[0][1];
    float largestY = v[0][1];
    
    float smallestZ = v[0][2];
    float largestZ = v[0][2];
    
    float avgX;
    float avgY;
    float avgZ;
    
    ofVec3f center;
    
    //Loop through v
    for (int i = 0; i < v.size(); i++)
    {
        //check and replace smallest values
        if(v[i][0] < smallestX) smallestX = v[i][0];
        if(v[i][1] < smallestY) smallestY = v[i][1];
        if(v[i][2] < smallestZ) smallestZ = v[i][2];
        
        //check and replace largest values
        if(v[i][0] > largestX) largestX = v[i][0];
        if(v[i][1] > largestY) largestY = v[i][1];
        if(v[i][2] > largestZ) largestZ = v[i][2];
    }
    
    //get averages from largest and smallest
    avgX = (smallestX + largestX)/2;
    avgY = (smallestY + largestY)/2;
    avgZ = (smallestZ + largestZ)/2;
    
    //create center vector
    center.set(avgX, avgY, avgZ);
    
    return center;
}

//--------------------------------------------------------------
void ofApp::clearFbo(ofFbo buffer)
{
    buffer.begin();
    ofClear(0,0,0,0);
    buffer.end();
}

//--------------------------------------------------------------
vector<string> ofApp::split(string str, char delimiter)
{
    vector<string> internal;
    stringstream ss(str);
    string tok;
    
    while(getline(ss, tok, delimiter))
    {
        internal.push_back(tok);
    }
    
    return internal;
}

//--------------------------------------------------------------
vector<float> ofApp::getFloatVector(vector<string> &sv)
{
    vector<float> floatVals;
    for(vector<string>::const_iterator itr = sv.begin(); itr != sv.end(); itr++)
    {
        //create string stream to check if float is valid
        stringstream sentry(*itr);
        float f; //stream data into float
        
        if(sentry >> f)
        {
            //if float is valid, push to float vector
            floatVals.push_back(f);
        }
    }
    
    return floatVals;
}

//--------------------------------------------------------------
string ofApp::removeSpaces(string str)
{
    int n = str.length();
    int i = 0, j = -1;
    bool spaceFound = false;
    
    //deal with leading spaces (if present)
    while (++j < n && str[j] == ' ');
    
    //read all characters of the original string
    while (j < n)
    {
        //if current characters are non-space characters
        if (str[j] != ' ')
        {
            //remove preceding spaces before dot, comma, & question mark
            if ((str[j] == '.' || str[j] == ',' || str[j] == '?') && str[i-1] == ' ')
            {
                str[i-1] = str[j++];
            }
            //copy current character at index i and increment both i and j
            else str[i++] = str[j++];
            
            //set space flag to false when any non-space character is found
            spaceFound = false;
        }
        // if current character is a space
        else if (str[j++] == ' ')
        {
            //If space is encounered for the first time after a word,
            //put one space in the outputand set space flag to true
            if(!spaceFound)
            {
                str[i++] = ' ';
                spaceFound = true;
            }
        }
    }
    //Remove trailing spaces
    if (i <= 1)
    {
        str.erase(str.begin() + i, str.end());
    }
    else
        str.erase(str.begin() + i - 1, str.end());
    
    return str;
    
}

//--------------------------------------------------------------
void ofApp::createNormals(ofMesh m)
{
    int vm = m.getNumVertices();
    for (int i = 0; i <= vm; i++)
    {
        //TODO: create normals for mesh
    }
}

//--------------------------------------------------------------
void ofApp::cycleSideChains(int indexNum)
{
    auto &verts = sideChainGroups[indexNum].getVertices();
    for(unsigned int i = 0; i < verts.size(); i++)
    {
        verts[i].x += ofRandom(-.1,.1), ofGetElapsedTimef();
        verts[i].y += ofRandom(-.1,.1), ofGetElapsedTimef();
        verts[i].z += ofRandom(-.1,.1), ofGetElapsedTimef();
    }
    sideChainGroups[indexNum].draw();
}

//--------------------------------------------------------------
void ofApp::drawBackbone()
{
    
    fbo.begin();
    ofClear(255,255,255,0);
    //Start Camera
    cam.begin();
    
    ofPushMatrix();
    //Move vertices to center
    ofTranslate(-10, 6, zVal);
    auto &verts = backbone.getVertices();
    for(unsigned int i = 0; i < verts.size(); i++)
    {
        verts[i].x += ofRandom(-.1,.1), ofGetElapsedTimef();
        verts[i].y += ofRandom(-.1,.1), ofGetElapsedTimef();
        verts[i].z += ofRandom(-.1,.1), ofGetElapsedTimef();
        //set vertex color
        //mesh.addColor(ofFloatColor(randRed,randGreen,randBlue,255));
        
        ofVec3f verta = backbone.getVertex(i);
        
        for (int j=i+1; j<verts.size(); ++j)
        {
            ofVec3f vertb = backbone.getVertex(j);
            //calculate distance between vertices
            float distance = verta.distance(vertb);
            //compare to threshold value (connectionDistance)
            //connectionDistance is set at top of draw()
            if (distance <= connectionDistance)
            {
                //set color for line
                ofSetLineWidth(2);
                //ofSetColor(233,253,255,150);
                ofSetColor(90,100,200,150);
                //ofSetColor(0,0,0,150);
                //draw line
                ofDrawLine(verts[i].x,verts[i].y,verts[i].z, vertb.x, vertb.y, vertb.z);
            }
        }
    }
    backbone.draw();
    
    //Set color for sidechains
    ofSetColor(250,10,0,255);
    
    //Cycle through the sidechain meshes
    cycleSideChains(mvCounter);
    
    ofPopMatrix();
    //Turn off camera
    cam.end();
    //End FBO
    fbo.end();
    
}

//--------------------------------------------------------------
void ofApp::drawMesh()
{
    //set drawtype
    //Drawtypes will be LINE, SLICE, POINTCLOUDS
    //int drawtype = type;
    //Take modulus of mesh size
    vertexCounter %= numVertices;
    fbo.begin();
    ofClear(255,255,255,0);
    //Start Camera
    cam.begin();
    
    ofPushMatrix();
    //Move vertices to center
    ofTranslate(-10, 6, zVal);
    auto &verts = mesh.getVertices();
    
    for(unsigned int i = 0; i < verts.size(); i++)
    {
        //Generate Noise
        verts[i].x += ofRandom(-.1,.1), ofGetElapsedTimef();
        verts[i].y += ofRandom(-.1,.1), ofGetElapsedTimef();
        verts[i].z += ofRandom(-.1,.1), ofGetElapsedTimef();
        
        //set vertex color
        //mesh.addColor(ofFloatColor(randRed,randGreen,randBlue,255));
    
        ofVec3f verta = mesh.getVertex(i);
        
        for (int j=i+1; j<verts.size(); ++j)
        {
            ofVec3f vertb = mesh.getVertex(j);
            //calculate distance between vertices
            float distance = verta.distance(vertb);
            //compare to threshold value (connectionDistance)
            //connectionDistance is set at top of draw()
            if (distance <= connectionDistance)
            {
                //set color for line
                ofSetLineWidth(1);
                //ofSetColor(233,253,255,150);
                ofSetColor(80,100,255,150);
                //draw line
                ofDrawLine(verts[i].x,verts[i].y,verts[i].z, vertb.x, vertb.y, vertb.z);
            }
        }
        
    }
    
    //Draw mesh
    mesh.draw();
    //ofSetColor(233,253,255,150);
    
    ofPopMatrix();
    
    //Turn off camera
    cam.end();
    //End FBO
    
    fbo.end();
    
}

//--------------------------------------------------------------
void ofApp::orbitCam()
{
    ofVec3f camPos = cam.getGlobalPosition();
    ofVec3f center = getCenter(coordinates);
    float theta = ofGetElapsedTimeMillis() * rotationSpeed;
    float x;
    float z;
    x = center.x + 100*cos(theta);
    z = center.z + 50*sin(theta);
    
    cam.setGlobalPosition(x, camPos.y, z);
    cam.setTarget(center);
}

//--------------------------------------------------------------
void ofApp::drawNoise(float noiseF)
{
    //
    //noiseFbo.begin();
    ofClear(0, 0, 0, 0);
    
    shaderNoise.begin();
    shaderNoise.setUniform1f("noiseAmt", noiseF);
    //TODO: set parameters for Bloom shader
    blurYFbo.draw(0, 0);
    shaderNoise.end();
    
    //noiseFbo.end();*/
}

//--------------------------------------------------------------
void ofApp::drawBloom(float intens)
{
    bloomFbo.begin();
    ofClear(255,255,255, 0);
    shaderBloom.begin();
    shaderBloom.setUniform1f("intensity", intens);
    fbo.draw(0,0);
    shaderBloom.end();
    bloomFbo.end();
}


//--------------------------------------------------------------
void ofApp::drawBlurX(float blr)
{
    blurXFbo.begin();
    ofClear(0,0,0,0);
    shaderBlurX.begin();
    shaderBlurX.setUniform1f("blurAmnt", blr);
    bloomFbo.draw(0, 0);
    shaderBlurX.end();
    blurXFbo.end();
}

//--------------------------------------------------------------
void ofApp::drawBlurY(float blr)
{
    blurYFbo.begin();
    ofClear(0, 0, 0, 0);
    shaderBlurY.begin();
    shaderBlurY.setUniform1f("blurAmnt", blr);
    blurXFbo.draw(0,0);
    shaderBlurY.end();
    blurYFbo.end();
}

//--------------------END USR FUNCTIONS-------------------------


//--------------------------------------------------------------
void ofApp::update()
{
    noiseF = ofRandom(-255, 255);
    
    //TODO: set messages with oscx.addFloatArg(xvalue); sender.sendMessage(oscx);
    mvCounter %= numSideChains;
    if (ofGetFrameNum() % 100 == 0)
    {
        tempo = int(ofRandom(5, 12));
    }
    //Handle OSC messaging
    //Create messages and set addresses
    ofxOscMessage oscNumSCAtoms;
    ofxOscMessage oscHydro;
    oscHydro.setAddress("/a");
    oscNumSCAtoms.setAddress("/b");
    
    
   //rotate model
    float fn = ofGetFrameNum();
    rotx = 0.03 * fn;
    roty = -0.02 * fn;
    rotz = 0.04 * fn;
    
    //Get current amino acid label and set to display string
    displayAminoAcid = aminoAcidLabels[mvCounter];
    std::map<string, string>::iterator it;
    it = hydrovals.find(displayAminoAcid);
    currHydro = it->second;
    
    //Add OSC arguments
    oscHydro.addFloatArg(stof(currHydro));
    oscNumSCAtoms.addIntArg(sideChainAtoms[mvCounter].size());
    //Send OSC messages
    sender.sendMessage(oscHydro, false);
    sender.sendMessage(oscNumSCAtoms, false);
    
    //Draw geometry into FBO
    //drawMesh(); //Draw to first FBO
    drawBackbone();
    drawBloom(intensity); //Draw to second FBO
    //drawNoise(); //Draw to third FBO
    drawBlurX(blurAmount);
    drawBlurY(blurAmount);
    
    //-------------------------------------------------------
    if(ofGetFrameNum() % tempo == 0)
    {
        displayVector.clear();
        for(std::vector<string>::const_iterator itr = sideChainAtoms[mvCounter].begin(); itr != sideChainAtoms[mvCounter].end(); ++itr)
        {
            //displayVector.resize(sideChainAtoms[mvCounter].size());
            displayVector.push_back(*itr);
        }
        
        //Increment counter
        mvCounter++;
        
    }

    if (ofGetFrameNum() % 150 == 0)
    {
        connectionDistance = ofRandom(2.65, 4.3);
        
        //connectionDistance = ofRandom(5, 7);
        //blurAmount = ofRandom(.5, 1.1);
    }
 
}

//--------------------------------------------------------------
void ofApp::draw()
{
    
    ofEnableDepthTest();
    orbitCam();
    //drawBlurY(blurAmount);
    drawNoise(noiseF);
    ofDisableDepthTest();
    //gui.draw();
    ofSetColor(100,100,100,255);
    ofDrawBitmapString("Glucagon-like peptide 1 receptor _ extracellular domain", 50, 60);
    ofDrawBitmapString("PDB _entry.id 5OTT", 50, 40);
    ofDrawBitmapString("amino acid", 50, 95);
    ofDrawBitmapString("sidechain molecule", 284, 95);
    ofDrawBitmapString("hydrophobicity", 152, 95);
    ofSetColor(250,10,0,255);
    ofDrawBitmapString(currHydro, 152, 115);
    ofSetColor(250,10,0,255);
    ofDrawBitmapString(displayAminoAcid, 50, 115);
    
    //counter for vector iterator
    int displayCounter=0;
   
   for(std::vector<string>::const_iterator itr = displayVector.begin(); itr != displayVector.end(); ++itr)
    {
        ofSetColor(250,10,0,255);
        ofDrawBitmapString((*itr), 284, 115+displayCounter*10);
        displayCounter++;
    }
    
    ////server.publishScreen();
    
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
