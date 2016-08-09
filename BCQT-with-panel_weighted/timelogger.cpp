#include "timelogger.h"

#include <iostream>
#include <iomanip>

TimeLogger* TimeLogger::timelogger = NULL;


TimeLogger* TimeLogger::Instance(){
  if(!timelogger){
    timelogger = new TimeLogger();
  }

  return timelogger;
}

void TimeLogger::push(QString label){
  timerinfos.top()->timer += timer->elapsed();
  timerinfos.push(new TimeInformation);
  if(label != "" && timerinfos.top() != NULL){
      Message * message = new Message(timerinfos.size(),label);
      messages.push_back(message);
      timerinfos.top()->messages.push_back(message);
      message->messagetype = MESSAGE;
  }

}

void TimeLogger::pop(QString label){
  TimeInformation* timeinfo = timerinfos.top();

  for(unsigned int i = 0; i < timeinfo->messages.size(); ++i){
     timeinfo->messages[i]->percent =  timeinfo->messages[i]->time*1.0/timeinfo->totaltime*100;
  }
  long time = timerinfos.top()->totaltime;

  if(label == "")
    markIt("Time to finish this level: ", time);
  else
    markIt(label, time);

  timerinfos.pop();

  if(timerinfos.size())
    timerinfos.top()->timer += time;
}

void TimeLogger::start(){
  if(timerinfos.size() == 0)
    timerinfos.push(new TimeInformation);
  timerinfos.top()->timer = 0;
  timer->start();
}


void TimeLogger::outputToScreen(){

  if(timerinfos.size())
    pop("The algorithm took: ");

  std::cout << std::setprecision(2) << std::fixed;
  if(messages.size()){
     for(int j = 1; j < messages[0]->level; ++j)
       std::cout<<"  ";

      if(messages[0]->messagetype & MESSAGE)
        std::cout<<messages[0]->message.toStdString();
      if(messages[0]->messagetype & TIME)
        std::cout<<convertTimeStr(messages[0]->time);
      if(messages[0]->messagetype & PERCENT)
        std::cout<<" %"<< messages[0]->percent;
      std::cout<<std::endl;

  }


  for(unsigned int i = 1; i < messages.size(); ++i){
    for(int j = 1; j < messages[i]->level; ++j)
      std::cout<<"  ";
    if(messages[i]->messagetype & MESSAGE)
      std::cout<<messages[i]->message.toStdString();
    if(messages[i]->messagetype & TIME)
      std::cout<<convertTimeStr(messages[i]->time);
    if(messages[i]->messagetype & PERCENT)
      std::cout<<" %"<< messages[i]->percent;
    std::cout<<std::endl;

  }


}

std::string TimeLogger::convertTimeStr(int msec){

  int hour, min, sec;

  if( msec > 1000*60*60){
     hour = msec/1000.0*60*60;
     min = (msec%(1000*60*60))/1000.0*60;
     sec = msec%(1000*60)/1000;
     return QString("%1h %2m %3s").arg(hour).arg(min).arg(sec).toStdString();
  }
  else if( msec > 1000*60){
      min = msec/1000.0*60;
      sec = msec%(1000*60)/1000;
      return QString("%1m %2s").arg(min).arg(sec).toStdString();
  }
  else if( msec > 1000){
      sec = msec/1000;

      return QString("%1s %2ms").arg(sec).arg(msec%1000).toStdString();
  }
  else if( msec > 0){
      return QString("%1ms").arg(msec).toStdString();
  }
  else{
      return QString("%1ms").arg(0).toStdString();
  }
}

std::string TimeLogger::elapsedTimeStr(){
 int timeelapsed = timer->elapsed();
 return convertTimeStr(timeelapsed);
}


float TimeLogger::elapsedTimerH(){
  return timer->elapsed()/1000.0/60/60;
}

float TimeLogger::elapsedTimerM(){
  return timer->elapsed()/1000.0/60;
}


float TimeLogger::elapsedTimeS(){
  return timer->elapsed()/1000.0;
}

int TimeLogger::elapsedTimeMs(){
  return timer->elapsed();
}

int TimeLogger::elapsedTimeNs(){
  return timer->nsecsElapsed();
}


void TimeLogger::markIt(QString str, int timems){
  Message * message = new Message(timerinfos.size(),str);
  messages.push_back(message);
  timerinfos.top()->messages.push_back(message);

  if(timems == -1){
    timems = timer->elapsed() + timerinfos.top()->timer;
    timerinfos.top()->totaltime += timems;
    start();
  }else{
   message->messagetype -= PERCENT;
  }
  message->time = timems;
}

TimeLogger::TimeLogger()
{
  timer = new QElapsedTimer();
}
