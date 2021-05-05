Events
======

Events are a way to react to something specific occuring, by adding a function known as a callback which will be called
each time the event is triggered. They are common in other languages such as C#.

*narupatools* uses events thouroughly to allow easy extensibility of code. Events allow you to insert code that runs for certain events such as after each dynamics step or at the start of each interaction. An event is created as so:

.. testcode::

   from narupatools.core.event import Event
   event = Event()

A callback is a function (normally a local function) that we wish to be called each time the event occurs. We add a callback to a event by either calling the `add_callback` function or by using the `+=` operator:

.. testcode::

   def callback(**kwargs):
       print("Hello World")

   event.add_callback(callback)

.. note::
   It is strongly advised that all callback functions have `**kwargs` as their final argument. This means they accept arbitrary arguments that will be ignored. If you don't do this, if the event gains a parameter in the future your callback will not work.

Adding the callback does not call the function. Also note we are passing in the name of the function (so `callback` instead of `callback()`).

Now, when we want to trigger the event, we call the `invoke()` method. This iterates through all the callbacks added to the event and calls them.

.. doctest::

   >>> event.invoke()
   Hello World
