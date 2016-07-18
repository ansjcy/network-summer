////
////  properties.hpp
////  BC
////
////  Created by Anakin on 16/7/14.
////  Copyright © 2016年 Anakin. All rights reserved.
////
//
//#ifndef properties_hpp
//#define properties_hpp
//
//#include <string>
//#include <stdio.h>
//#include <stdlib.h>
//#include <vector>
//#include <map>
//using namespace std;
//
//typedef enum {
//    PROPERTY_UNKNOWN,
//    PROPERTY_BOOL,
//    PROPERTY_INT,
//    PROPERTY_FLOAT,
//    PROPERTY_SIGNED,
//    PROPERTY_VEC2,
//    PROPERTY_VEC3,
//    PROPERTY_STRING,
//    PROPERTY_SELECTION,
//    PROPERTY_MULTISELECTION
//} PropType;
//
//class Property {
//protected:
//    PropType type;
//public:
//    Property() { type = PROPERTY_UNKNOWN; }
//    Property(PropType t) { type = t; }
//    PropType getType() { return type;}
//    virtual ~Property() { } //printf("Property type %d.Destroy\n", type);}
//    virtual string toString() { return ""; }
//    virtual void parseString(string value) { };
//};
//
//class StringProperty:public Property {
//public:
//    string value;
//    StringProperty(const char* s): Property(PROPERTY_STRING) {
//        this->value = string(s);
//    }
//    StringProperty(string  s): Property(PROPERTY_STRING) {
//        this->value = s;
//    }
//    virtual string toString() { return value; }
//    virtual void parseString(string s) { this->value = s; }
//};
//
//class ListProperty:public StringProperty {
//public:
//    string model;
//    ListProperty(const char* s, const char *model = 0): StringProperty(s) {
//        this->model = model? string(model): "none";
//        this->type = PROPERTY_SELECTION;
//    }
//    ListProperty(string s, string model = "none"): StringProperty(s) {
//        this->model = model;
//        this->type = PROPERTY_SELECTION;
//    }
//    virtual string toString() { return value; }
//    virtual void parseString(string s) { this->value = s; }
//};
//
//class MultiListProperty:public StringProperty {
//public:
//    string model;
//    MultiListProperty(const char* s, const char *model = 0): StringProperty(s) {
//        this->model = model? string(model): "none";
//        this->type = PROPERTY_MULTISELECTION;
//    }
//    MultiListProperty(string s, string model = "none"): StringProperty(s) {
//        this->model = model;
//        this->type = PROPERTY_MULTISELECTION;
//    }
//    virtual string toString() { return value; }
//    virtual void parseString(string s) { this->value = s; }
//};
//
//class BoolProperty:public Property {
//public:
//    bool value;
//    BoolProperty(bool v): Property(PROPERTY_BOOL) {
//        this->value = v;
//    }
//    BoolProperty(BoolProperty &p): Property(PROPERTY_BOOL) {
//        this->value = p.value;
//    }
//    virtual string toString() {
//        return string(value? "true":"false");
//    }
//    virtual void parseString(string s) { this->value = (s=="true"); }
//};
//
//class IntProperty:public Property {
//public:
//    int value;
//    IntProperty(int v): Property(PROPERTY_INT) {
//        this->value = v;
//    }
//    IntProperty(IntProperty &p): Property(PROPERTY_INT) {
//        this->value = p.value;
//    }
//    virtual string toString() {
//        char str[200];
//        sprintf(str, "%d", value);
//        return string(str);
//    }
//    virtual void parseString(string s) { this->value = atoi(s.c_str()); }
//};
//
//class FloatProperty:public Property {
//public:
//    float value;
//    FloatProperty(float v): Property(PROPERTY_FLOAT) {
//        this->value = v;
//    }
//    FloatProperty(FloatProperty &p): Property(PROPERTY_FLOAT) {
//        this->value = p.value;
//    }
//    virtual string toString() {
//        char str[200];
//        sprintf(str, "%f", value);
//        return string(str);
//    }
//    virtual void parseString(string s) { this->value = atof(s.c_str()); }
//};
//
//class SignedFloatProperty:public FloatProperty {
//public:
//    SignedFloatProperty(float v): FloatProperty(v) {
//        type = PROPERTY_SIGNED;
//    }
//    SignedFloatProperty(SignedFloatProperty &p): FloatProperty(p.value) {
//        type = PROPERTY_SIGNED;
//    }
//};
//
//class Vec2Property: public Property {
//public:
//    float vector[2];
//    Vec2Property(float x, float y): Property(PROPERTY_VEC2) {
//        vector[0] = x;
//        vector[1] = y;
//    }
//    Vec2Property(float * v): Property(PROPERTY_VEC2) {
//        vector[0] = v[0];
//        vector[1] = v[1];
//    }
//    float operator[](int i) {
//        return vector[i];
//    }
//    virtual string toString() {
//        char str[200];
//        sprintf(str, "(%f,%f)", vector[0], vector[1]);
//        return string(str);
//    }
//    virtual void parseString(string s) {
//        sscanf(s.c_str(), "(%f,%f)", &(vector[0]), &(vector[1]));
//    }
//};
//
//class Vec3Property: public Property {
//public:
//    float vector[3];
//    Vec3Property(float x, float y, float z): Property(PROPERTY_VEC3) {
//        vector[0] = x;
//        vector[1] = y;
//        vector[2] = z;
//    }
//    Vec3Property(float *v): Property(PROPERTY_VEC3) {
//        vector[0] = v[0];
//        vector[1] = v[1];
//        vector[2] = v[2];
//    }
//    float operator[](int i) {
//        return vector[i];
//    }
//    virtual string toString() { 
//        char str[200];
//        sprintf(str, "(%f,%f,%f)", vector[0], vector[1], vector[2]);
//        return string(str); 
//    }
//    virtual void parseString(string s) { 
//        sscanf(s.c_str(), "(%f,%f,%f)", &(vector[0]), &(vector[1]), &(vector[2]));
//    }
//};
//
//
//using namespace std;
//
//class PropertyBag{
//protected:
//    map<string,Property*> properties;
//public:
//    PropertyBag();
//    ~PropertyBag();
//    vector<string> getProperties();
//    bool hasProperty(string name);
//    void setProperty(string name, Property *p);
//    void setProperty(string name, string value);
//    Property* getProperty(string name);
//    
//    void setPropertyString(string name, string s);
//    void setPropertyString(string name, const char * s);
//    
//    void setPropertyStringList(string name, string s, string model = "none", bool multi=false);
//    void setPropertyStringList(string name, const char * s, const char * model = 0, bool multi=false);
//    
//    void setPropertyFloat(string name, float v);
//    void setPropertySigned(string name, float v);
//    void setPropertyInt(string name, int v);
//    void setPropertyBool(string name, bool v);
//    void setPropertyVec2(string name, float x, float y);
//    void setPropertyVec2v(string name, float *v);
//    void setPropertyVec3(string name, float x, float y, float z);
//    void setPropertyVec3v(string name, float *v);
//    
//    string getPropertyString(string name);
//    int getPropertyInt(string name, int defaultValue = 0);
//    bool getPropertyBool(string name, bool defaultValue = false);
//    float getPropertyFloat(string name, float defaultValue = 0.0);
//    void getPropertyVec2(string name, float *v);
//    void getPropertyVec3(string name, float *v);
//    
//    int getPropertyType(string name);
//    string getPropertyListModel(string name);
//    
//};
//#endif /* properties_hpp */
