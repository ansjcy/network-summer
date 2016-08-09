#ifndef TIMELOGGER_H
#define TIMELOGGER_H

#include <QElapsedTimer>
#include <QString>
#include <QStringList>
#include <QList>
#include <QStack>

enum message_t{
  FULL_MESSAGE = 7,
  TIME = 4,
  MESSAGE = 2,
  PERCENT = 1
};


struct Message{

   Message(int level, QString message){
    this->level = level;
    this->message = message;
    this->percent = 0;
    this->time = 0;
    this->messagetype = FULL_MESSAGE;
   }

   int level;
   long time;
   float percent;
   int messagetype;
   QString message;
};

struct TimeInformation{


  TimeInformation(){
    timer = 0;
    totaltime = 0;
  }

  long timer;
  std::vector<Message*> messages;
  long totaltime;
};

class TimeLogger
{
public:
  static TimeLogger* Instance();

  void push(QString label = "");
  void pop(QString label = "");

  void start();
  void markIt(QString str, int timems = -1);

  void outputToScreen();
  void outputToFile(QString *file);

  std::string elapsedTimeStr();

  float elapsedTimerH();
  float elapsedTimerM();
  float elapsedTimeS();
  int elapsedTimeMs();
  int elapsedTimeNs();


private:
  TimeLogger();
  std::string convertTimeStr(int msec);

  static TimeLogger* timelogger;

  QElapsedTimer *timer;

  QStack<TimeInformation *> timerinfos;
  std::vector<Message*> messages;

};

#endif // TIMELOGGER_H
