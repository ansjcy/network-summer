////
////  properties.cpp
////  BC
////
////  Created by Anakin on 16/7/14.
////  Copyright © 2016年 Anakin. All rights reserved.
////
//
//#include "properties.hpp"
//
//PropertyBag::~PropertyBag() {
//    // print("PropertyBag.Destroy\n");
//    
//    for(map<string,Property*>::const_iterator it = properties.begin(); it !=properties.end(); ++it) {
//        delete (it->second);
//    }
//    
//    properties.clear();
//}
//
//vector<string> PropertyBag::getProperties() {
//    vector<string> vec;
//    //printf("PropertyBag.Destroy\n");
//    for(map<string,Property*>::const_iterator it = properties.begin(); it !=properties.end(); ++it) {
//        vec.push_back((it->first));
//    }
//    return vec;
//}
//bool PropertyBag::hasProperty(string name) {
//    return (properties.find(name)!=properties.end());
//}
//
//void PropertyBag::setProperty(string name, Property *p) {
//    //print("Adding property %s %s\n", name.c_str(), p->toString().c_str());
//    properties[name] = p;
//}
//
//void PropertyBag::setProperty(string name, string value) {
//    Property *prop = getProperty(name);
//    if(prop) {
//        prop->parseString(value);
//    }
//}
//
//Property* PropertyBag::getProperty(string name) {
//    if(!hasProperty(name)) return 0;
//    return properties[name];
//}
//
//void PropertyBag::setPropertyString(string name, string s) {
//    properties[name] = new StringProperty(s);
//}
//void PropertyBag::setPropertyString(string name, const char* s) {
//    properties[name] = new StringProperty(s);
//}
//
//void PropertyBag::setPropertyStringList(string name, string s, string model, bool multi) {
//    if(multi) {
//        properties[name] = new MultiListProperty(s, model);
//    } else {
//        properties[name] = new ListProperty(s, model);
//    }
//}
//void PropertyBag::setPropertyStringList(string name, const char* s, const char *model, bool multi) {
//    if(multi) {
//        properties[name] = new MultiListProperty(s, model);
//    } else {
//        properties[name] = new ListProperty(s, model);
//    }
//}
//
//
//void PropertyBag::setPropertyFloat(string name, float v) {
//    properties[name] = new FloatProperty(v);
//    //print("Added property value %s=%f\n", name.c_str(), v);
//}
//
//void PropertyBag::setPropertySigned(string name, float v) {
//    //print("Added property signed value %s=%f\n", name.c_str(), v);
//    properties[name] = new SignedFloatProperty(v);
//}
//
//void PropertyBag::setPropertyInt(string name, int v) {
//    properties[name] = new IntProperty(v);
//}
//
//void PropertyBag::setPropertyBool(string name, bool v) {
//    properties[name] = new BoolProperty(v);
//}
//
//
//void PropertyBag::setPropertyVec2(string name, float x, float y) {
//    properties[name] = new Vec2Property(x,y);
//}
//
//void PropertyBag::setPropertyVec2v(string name, float *v) {
//    properties[name] = new Vec2Property(v);
//}
//
//void PropertyBag::setPropertyVec3(string name, float x, float y, float z) {
//    properties[name] = new Vec2Property(x,y);
//}
//
//void PropertyBag::setPropertyVec3v(string name, float *v) {
//    properties[name] = new Vec3Property(v);
//}
//
//float PropertyBag::getPropertyFloat(string name, float defaultValue) {
//    if(properties.find(name)==properties.end()) {
//        //error("Could not find property %s\n", name.c_str());
//        return defaultValue;
//    }
//    FloatProperty* p = (FloatProperty*) (properties[name]);
//    return p->value;
//}
//
//int PropertyBag::getPropertyInt(string name, int defaultValue) {
//    if(properties.find(name)==properties.end()) {
//        return defaultValue;
//    }
//    IntProperty* p = (IntProperty*) (properties[name]);
//    return p->value;
//}
//
//bool PropertyBag::getPropertyBool(string name, bool defaultValue) {
//    if(properties.find(name)==properties.end()) {
//        return defaultValue;
//    }
//    BoolProperty* p = (BoolProperty*) (properties[name]);
//    return p->value;
//}
//
//void PropertyBag::getPropertyVec2(string name, float *v) {
//    Vec2Property* p = (Vec2Property*) properties[name];
//    v[0] = p->vector[0];
//    v[1] = p->vector[1];
//}
//
//void PropertyBag::getPropertyVec3(string name, float *v) {
//    Vec3Property* p = (Vec3Property*) (properties[name]);
//    v[0] = p->vector[0];
//    v[1] = p->vector[1];
//    v[2] = p->vector[2];
//}
//
//string PropertyBag::getPropertyString(string name) {
//    if(properties.find(name)==properties.end()) {
//        return "";
//    }
//    return ((StringProperty*) properties[name])->value;
//}
//
//
//int PropertyBag::getPropertyType(string name) {
//    if(properties.find(name)==properties.end()) {
//        return PROPERTY_UNKNOWN;
//    }
//    return properties[name]->getType();
//}
//
//string PropertyBag::getPropertyListModel(string name) {
//    Property *p = properties[name];
//    if(!p) return "";
//    if(p->getType()==PROPERTY_SELECTION) {
//        return ((ListProperty*) p)->model;
//    }
//    if(p->getType()==PROPERTY_MULTISELECTION) {
//        return ((MultiListProperty*) p)->model;
//    }
//    return "";
//}
