"""Explicit observable FLINT-call accounting, with no FLINT import or numerical probe."""
import time,traceback
CATEGORIES=('setup','construction','arithmetic','query')
def empty_counts():return {'attempted':0,'completed':0,'raised_exception':0}
def new_ledger():return {'totals':empty_counts(),'by_category':{k:empty_counts()for k in CATEGORIES},'budgeted_setup_construction_arithmetic_attempts':0,'budget_rejections':0,'inflight_api':None,'last_raised_exception':None,'observed_stage':None,'exactness':'observed_in_this_live_process; process death makes uncheckpointed totals UNKNOWN'}
class Ledger:
 def __init__(self,global_state,caps,stage_label):
  self.global_state=global_state;self.caps=caps;self.stage=new_ledger();global_state['observed_stage']={'label':stage_label,'ledger':self.stage}
 def invoke(self,category,label,fn):
  if category not in CATEGORIES:raise ArithmeticError('UNDECLARED_API_CATEGORY')
  g=self.global_state;s=self.stage
  if category!='query'and(s['budgeted_setup_construction_arithmetic_attempts']>=self.caps['native_arithmetic_per_stage']or g['budgeted_setup_construction_arithmetic_attempts']>=self.caps['native_arithmetic_total']):
   s['budget_rejections']+=1;g['budget_rejections']+=1;raise ArithmeticError('SETUP_CONSTRUCTION_ARITHMETIC_ATTEMPT_CAP')
  event={'global_attempt_id':g['totals']['attempted']+1,'category':category,'label':label,'state':'INFLIGHT_UNKNOWN','monotonic_start':time.monotonic()}
  # Mark uncertainty FIRST. An interrupt anywhere in subsequent bookkeeping remains UNKNOWN.
  g['inflight_api']=dict(event);s['inflight_api']=dict(event)
  if category!='query':
   s['budgeted_setup_construction_arithmetic_attempts']+=1;g['budgeted_setup_construction_arithmetic_attempts']+=1
  for state in (s,g):
   state['totals']['attempted']+=1;state['by_category'][category]['attempted']+=1
  try:result=fn()
  except BaseException:
   failed={**event,'state':'INVOCATION_RAISED_EXCEPTION','exception':traceback.format_exc()}
   for state in (s,g):
    state['totals']['raised_exception']+=1;state['by_category'][category]['raised_exception']+=1;state['last_raised_exception']=failed;state['inflight_api']=None
   raise
  for state in (s,g):
   state['totals']['completed']+=1;state['by_category'][category]['completed']+=1;state['inflight_api']=None
  return result
 def check_quiescent(self):
  for state in (self.stage,self.global_state):
   if state['inflight_api']is not None or state['totals']['attempted']!=state['totals']['completed']+state['totals']['raised_exception']:raise ArithmeticError('API_LEDGER_QUIESCENT_IDENTITY')
