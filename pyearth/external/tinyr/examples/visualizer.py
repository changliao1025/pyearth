#
#    tinyr - a 2D-RTree implementation in Cython
#    Copyright (C) 2011  Matthias Simon
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



'''
Interactive visualizer of tinyr RTree.

Displays inserted leafs of RTree and the surrounding index rectangle of leafs
(RTree internal). 

Usage for now:
    - space bar inserts a random rectangle
    - mouse wheel zooms in/out
    - left-click drag drags the display
    - right-click drag draws a rectangle, when done: overlapping leaf nodes
      are removed


Most code for visualization and user action is taken from xdot (Interactive 
viewer for Graphviz dot files, http://code.google.com/p/jrfonseca/wiki/XDot )

'''

import math
import time
import colorsys

import gtk
import gobject

from util import tinyr, node_gen

class Rectangle(object):
    def __init__(self, rect, color=(0.0, 0.0, 0.0, 1.0), linewidth=2, fillcolor=None):
        self.rect = rect
        self.color = color
        self.linewidth = linewidth
        self.fillcolor = fillcolor
    
    def draw(self, cr):
        
        cr.save()
        cr.set_source_rgba(*self.color)
        cr.set_line_width(self.linewidth)

        cr.rectangle(*self.rect)
        cr.stroke()
        
        if self.fillcolor:
            cr.set_source_rgba(*pen.fillcolor)
            cr.fill_preserve()
            cr.fill()

        cr.restore()
        
class ItemRectangle(Rectangle):
    def __init__(self, rect, level, levels_cnt):
        self.rect = rect
        
        if level == levels_cnt+1:
            color = (0.0, 0.0, 0.0, 1.0)
        else:
            levels_cnt = max(1, levels_cnt)
            color = list(colorsys.hsv_to_rgb(float(level+1)/levels_cnt, 1.0, 1.0))
            color.append(1.0)
        
        Rectangle.__init__(self, rect, color, (level+1) * 2)
        self.level = level
    
    def draw(self, cr):
        
        cr.save()
        cr.set_source_rgba(*self.color)
        cr.set_line_width(self.linewidth)

        cr.rectangle(*self.rect)
        cr.stroke()
        cr.restore()


class Animation(object):
    step = 0.03 # seconds
    def __init__(self, rtree_widget):
        self.rtree_widget = rtree_widget
        self.timeout_id = None
    def start(self):
        self.timeout_id = gobject.timeout_add(int(self.step * 1000), self.tick)
    def stop(self):
        self.rtree_widget.animation = NoAnimation(self.rtree_widget)
        if self.timeout_id is not None:
            gobject.source_remove(self.timeout_id)
            self.timeout_id = None
    def tick(self):
        self.stop()

class NoAnimation(Animation):
    def start(self):
        pass
    def stop(self):
        pass

class LinearAnimation(Animation):
    duration = 0.6
    def start(self):
        self.started = time.time()
        Animation.start(self)
    def tick(self):
        t = (time.time() - self.started) / self.duration
        self.animate(max(0, min(t, 1)))
        return (t < 1)
    def animate(self, t):
        pass

class MoveToAnimation(LinearAnimation):
    def __init__(self, rtree_widget, target_x, target_y):
        Animation.__init__(self, rtree_widget)
        self.source_x = rtree_widget.x
        self.source_y = rtree_widget.y
        self.target_x = target_x
        self.target_y = target_y
    def animate(self, t):
        sx, sy = self.source_x, self.source_y
        tx, ty = self.target_x, self.target_y
        self.rtree_widget.x = tx * t + sx * (1-t)
        self.rtree_widget.y = ty * t + sy * (1-t)
        self.rtree_widget.queue_draw()


class ZoomToAnimation(MoveToAnimation):
    def __init__(self, rtree_widget, target_x, target_y):
        MoveToAnimation.__init__(self, rtree_widget, target_x, target_y)
        self.source_zoom = rtree_widget.zoom_ratio
        self.target_zoom = self.source_zoom
        self.extra_zoom = 0
        middle_zoom = 0.5 * (self.source_zoom + self.target_zoom)
        distance = math.hypot(self.source_x - self.target_x,
                              self.source_y - self.target_y)
        rect = self.rtree_widget.get_allocation()
        visible = min(rect.width, rect.height) / self.rtree_widget.zoom_ratio
        visible *= 0.9
        if distance > 0:
            desired_middle_zoom = visible / distance
            self.extra_zoom = min(0, 4 * (desired_middle_zoom - middle_zoom))

    def animate(self, t):
        a, b, c = self.source_zoom, self.extra_zoom, self.target_zoom
        self.rtree_widget.zoom_ratio = c*t + b*t*(1-t) + a*(1-t)
        self.rtree_widget.zoom_to_fit_on_resize = False
        MoveToAnimation.animate(self, t)


class DragAction(object):
    def __init__(self, rtree_widget):
        self.rtree_widget = rtree_widget
    def on_button_press(self, event):
        self.startmousex = self.prevmousex = event.x
        self.startmousey = self.prevmousey = event.y
        self.start()
    def on_motion_notify(self, event):
        deltax = self.prevmousex - event.x
        deltay = self.prevmousey - event.y
        self.drag(deltax, deltay)
        self.prevmousex = event.x
        self.prevmousey = event.y
    def on_button_release(self, event):
        self.stopmousex = event.x
        self.stopmousey = event.y
        self.stop()
    def draw(self, cr):
        pass

    def start(self):
        pass

    def drag(self, deltax, deltay):
        pass

    def stop(self):
        pass

    def abort(self):
        pass


class NullAction(DragAction):

    def on_motion_notify(self, event):
        self.rtree_widget.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.ARROW))


class PanAction(DragAction):

    def start(self):
        self.rtree_widget.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.FLEUR))

    def drag(self, deltax, deltay):
        self.rtree_widget.x += deltax / self.rtree_widget.zoom_ratio
        self.rtree_widget.y += deltay / self.rtree_widget.zoom_ratio
        self.rtree_widget.queue_draw()

    def stop(self):
        self.rtree_widget.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.ARROW))

    abort = stop

    abort = stop

class RectAction(DragAction):

    def drag(self, deltax, deltay):
        self.rtree_widget.queue_draw()

    def draw(self, cr):
        cr.save()
        cr.set_source_rgba(0.0, 0.0, 1.0, 1.0)
        cr.set_line_width(2)

        x1, y1 = self.startmousex, self.startmousey
        x2, y2 = self.prevmousex, self.prevmousey

        cr.rectangle(x1, y1, x2-x1, y2-y1)

        cr.set_source_rgba(0.0, 0.0, 1.0, 0.2)
        cr.fill_preserve()
        cr.fill()

        cr.restore()

    def stop(self):
        x1, y1 = self.rtree_widget.window2rtree(self.startmousex,
                                              self.startmousey)
        x2, y2 = self.rtree_widget.window2rtree(self.stopmousex,
                                              self.stopmousey)
        
        identds = self.rtree_widget.rtree.search((x1, y1, x2, y2))
        for i in identds:
            self.rtree_widget.rtree.remove(i, i)
        
        self.rtree_widget.update_rectangles()
        self.rtree_widget.queue_draw()

class ZoomAction(DragAction):

    def drag(self, deltax, deltay):
        self.rtree_widget.zoom_ratio *= 1.005 ** (deltax + deltay)
        self.rtree_widget.zoom_to_fit_on_resize = False
        self.rtree_widget.queue_draw()

    def stop(self):
        self.rtree_widget.queue_draw()


class ZoomAreaAction(DragAction):

    def drag(self, deltax, deltay):
        self.rtree_widget.queue_draw()

    def draw(self, cr):
        cr.save()
        cr.set_source_rgba(.5, .5, 1.0, 0.25)
        cr.rectangle(self.startmousex, self.startmousey,
                     self.prevmousex - self.startmousex,
                     self.prevmousey - self.startmousey)
        cr.fill()
        cr.set_source_rgba(.5, .5, 1.0, 1.0)
        cr.set_line_width(1)
        cr.rectangle(self.startmousex - .5, self.startmousey - .5,
                     self.prevmousex - self.startmousex + 1,
                     self.prevmousey - self.startmousey + 1)
        cr.stroke()
        cr.restore()

    def stop(self):
        x1, y1 = self.rtree_widget.window2rtree(self.startmousex,
                                              self.startmousey)
        x2, y2 = self.rtree_widget.window2rtree(self.stopmousex,
                                              self.stopmousey)
        self.rtree_widget.zoom_to_area(x1, y1, x2, y2)

    def abort(self):
        self.rtree_widget.queue_draw()


class RTreeWidget(gtk.DrawingArea):
    __gsignals__ = {
        'expose-event': 'override',
        'clicked' : (gobject.SIGNAL_RUN_LAST, gobject.TYPE_NONE, (gobject.TYPE_STRING, gtk.gdk.Event))
    }

    def __init__(self):
        gtk.DrawingArea.__init__(self)

        self.rtree = None
        self.rectangles = []

        self.set_flags(gtk.CAN_FOCUS)

        self.add_events(gtk.gdk.BUTTON_PRESS_MASK | gtk.gdk.BUTTON_RELEASE_MASK)
        self.connect("button-press-event", self.on_area_button_press)
        self.connect("button-release-event", self.on_area_button_release)
        self.add_events(gtk.gdk.POINTER_MOTION_MASK | gtk.gdk.POINTER_MOTION_HINT_MASK | gtk.gdk.BUTTON_RELEASE_MASK)
        self.connect("motion-notify-event", self.on_area_motion_notify)
        self.connect("scroll-event", self.on_area_scroll_event)
        self.connect("size-allocate", self.on_area_size_allocate)

        self.x, self.y = 0.0, 0.0
        self.zoom_ratio = 1.0
        self.zoom_to_fit_on_resize = False
        self.animation = NoAnimation(self)
        self.drag_action = NullAction(self)
        self.presstime = None
        self.highlight = None

    def do_expose_event(self, event):
        cr = self.window.cairo_create()

        # set a clip region for the expose event
        cr.rectangle(
            event.area.x, event.area.y,
            event.area.width, event.area.height)
        cr.clip()

        cr.set_source_rgba(1.0, 1.0, 1.0, 1.0)
        cr.paint()

        cr.save()
        rect = self.get_allocation()
        cr.translate(0.5*rect.width, 0.5*rect.height)
        cr.scale(self.zoom_ratio, self.zoom_ratio)
        cr.translate(-self.x, -self.y)

        for rect in self.rectangles:
            rect.draw(cr)

        cr.restore()
        self.drag_action.draw(cr)

        return False

    def zoom_image(self, zoom_ratio, center=False, pos=None):
        if center:
            self.x = self.retreeinfo.width/2
            self.y = self.retreeinfo.height/2
        elif pos is not None:
            rect = self.get_allocation()
            x, y = pos
            x -= 0.5*rect.width
            y -= 0.5*rect.height
            self.x += x / self.zoom_ratio - x / zoom_ratio
            self.y += y / self.zoom_ratio - y / zoom_ratio
        self.zoom_ratio = zoom_ratio
        self.zoom_to_fit_on_resize = False
        self.queue_draw()

    def zoom_to_area(self, x1, y1, x2, y2):
        rect = self.get_allocation()
        width = abs(x1 - x2)
        height = abs(y1 - y2)
        self.zoom_ratio = min(
            float(rect.width)/float(width),
            float(rect.height)/float(height)
        )
        self.zoom_to_fit_on_resize = False
        self.x = (x1 + x2) / 2
        self.y = (y1 + y2) / 2
        self.queue_draw()

    def zoom_to_fit(self):
        rect = self.get_allocation()
        rect.x += self.ZOOM_TO_FIT_MARGIN
        rect.y += self.ZOOM_TO_FIT_MARGIN
        rect.width -= 2 * self.ZOOM_TO_FIT_MARGIN
        rect.height -= 2 * self.ZOOM_TO_FIT_MARGIN
        zoom_ratio = min(
            float(rect.width)/float(self.retreeinfo.width),
            float(rect.height)/float(self.retreeinfo.height)
        )
        self.zoom_image(zoom_ratio, center=True)
        self.zoom_to_fit_on_resize = True

    ZOOM_INCREMENT = 1.25
    ZOOM_TO_FIT_MARGIN = 12

    def on_zoom_in(self, action):
        self.zoom_image(self.zoom_ratio * self.ZOOM_INCREMENT)

    def on_zoom_out(self, action):
        self.zoom_image(self.zoom_ratio / self.ZOOM_INCREMENT)

    def on_zoom_fit(self, action):
        self.zoom_to_fit()

    def on_zoom_100(self, action):
        self.zoom_image(1.0)

    POS_INCREMENT = 100

    def get_drag_action(self, event):
        state = event.state
        if event.button in (1, 2): # left or middle button
            if state & gtk.gdk.CONTROL_MASK:
                return ZoomAction
            elif state & gtk.gdk.SHIFT_MASK:
                return ZoomAreaAction
            else:
                return PanAction
        elif event.button == 3:
            return RectAction
        return NullAction

    def on_area_button_press(self, area, event):
        self.animation.stop()
        self.drag_action.abort()
        action_type = self.get_drag_action(event)
        self.drag_action = action_type(self)
        self.drag_action.on_button_press(event)
        self.presstime = time.time()
        self.pressx = event.x
        self.pressy = event.y
        return False

    def is_click(self, event, click_fuzz=4, click_timeout=1.0):
        assert event.type == gtk.gdk.BUTTON_RELEASE
        if self.presstime is None:
            return False
        deltax = self.pressx - event.x
        deltay = self.pressy - event.y
        return (time.time() < self.presstime + click_timeout
                and math.hypot(deltax, deltay) < click_fuzz)

    def on_area_button_release(self, area, event):
        self.drag_action.on_button_release(event)
        self.drag_action = NullAction(self)
        if event.button == 1 and self.is_click(event):
            return True
        if event.button == 1 or event.button == 2:
            return True
        return False

    def on_area_scroll_event(self, area, event):
        if event.direction == gtk.gdk.SCROLL_UP:
            self.zoom_image(self.zoom_ratio * self.ZOOM_INCREMENT,
                            pos=(event.x, event.y))
            return True
        if event.direction == gtk.gdk.SCROLL_DOWN:
            self.zoom_image(self.zoom_ratio / self.ZOOM_INCREMENT,
                            pos=(event.x, event.y))
            return True
        return False

    def on_area_motion_notify(self, area, event):
        self.drag_action.on_motion_notify(event)
        return True

    def on_area_size_allocate(self, area, allocation):
        if self.zoom_to_fit_on_resize:
            self.zoom_to_fit()

    def animate_to(self, x, y):
        self.animation = ZoomToAnimation(self, x, y)
        self.animation.start()

    def window2rtree(self, x, y):
        rect = self.get_allocation()
        x -= 0.5*rect.width
        y -= 0.5*rect.height
        x /= self.zoom_ratio
        y /= self.zoom_ratio
        x += self.x
        y += self.y
        return x, y

    def set_rtree(self, rtree):
        self.rtree = rtree
        self.retreeinfo = rtree.get_info()
        self.update_rectangles()
    
    def update_rectangles(self):
        self.rectangles = []
        for level, rect in self.retreeinfo.iter_rectangles():
            rect = [ rect[0], rect[1], rect[2]-rect[0], rect[3]-rect[1] ]
            self.rectangles.append(ItemRectangle(rect, level, self.retreeinfo.levels))
        self.rectangles.sort(key=lambda x: x.level, reverse=True)
        self.queue_draw()

class RTreeWindow(gtk.Window):

    def __init__(self, rtree, size=100, shift=1200):
        gtk.Window.__init__(self)
        self.connect('key_press_event', self.on_key)
        
        self.rectangles = []

        window = self

        window.set_title('RTree Visualizer')
        window.set_default_size(800, 600)
        vbox = gtk.VBox()
        window.add(vbox)

        self.widget = RTreeWidget()
        uimanager = self.uimanager = gtk.UIManager()
        accelgroup = uimanager.get_accel_group()
        window.add_accel_group(accelgroup)

        actiongroup = gtk.ActionGroup('Actions')
        self.actiongroup = actiongroup

        actiongroup.add_actions((
            ('ZoomIn', gtk.STOCK_ZOOM_IN, None, None, None, self.widget.on_zoom_in),
            ('ZoomOut', gtk.STOCK_ZOOM_OUT, None, None, None, self.widget.on_zoom_out),
            ('ZoomFit', gtk.STOCK_ZOOM_FIT, None, None, None, self.widget.on_zoom_fit),
            ('Zoom100', gtk.STOCK_ZOOM_100, None, None, None, self.widget.on_zoom_100),
        ))

        uimanager.insert_action_group(actiongroup, 0)

        vbox.pack_start(self.widget)
        self.set_focus(self.widget)
        self.show_all()
        
        self.rtree = rtree
        self.nodegen = node_gen(99999, size=size, shift=shift)
        
        self.widget.zoom_to_area(0, 0, shift+size, shift+size)
        
        self.on_key(None, None)
    
    def on_key(self, win, ev):
        if ev and ev.keyval != 32:
            return True
        
        ident, coords = self.nodegen.next()
        self.rtree.insert(coords, coords)
        self.widget.set_rtree(self.rtree)
        self.widget.queue_draw()
    

def main():

    rt = tinyr.RTree()
    win = RTreeWindow(rt)
    
    win.connect('destroy', gtk.main_quit)
    gtk.main()


if __name__ == '__main__':
    main()




