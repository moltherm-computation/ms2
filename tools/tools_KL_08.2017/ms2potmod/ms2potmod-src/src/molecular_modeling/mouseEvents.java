/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package molecular_modeling;

/**
 *
 * @author Syed ahsan ali
 * 
 * This class is used to register mouse events inside scenegraph
 * - Get the specified location where user clicked..
 * - call function to add atom on scenegraph when clicked on appropriate atom...
 */
import com.sun.j3d.utils.behaviors.picking.Intersect;
import com.sun.j3d.utils.picking.*;

import com.sun.j3d.utils.universe.SimpleUniverse;

import com.sun.j3d.utils.geometry.*;

import javax.media.j3d.*;

import javax.vecmath.*;

import java.awt.event.*;

import java.awt.*;
import javax.swing.JOptionPane;

/**
 *
 * @author Syed Ahsan Ali
 */
public class mouseEvents extends MouseAdapter {

private PickCanvas pickCanvas;
private CreateGraphics Cmain;
Vector3d old = new Vector3d();
Shape3D s = null;
Primitive p= null;
Canvas3D canvas3d;

    /**
     * Add mouse listners.
     * 
     * @param group Reference of branch group
     * @param canvas Reference of canvas
     * @param parent Reference of CreateGraphics class
     * 
     * @see CreateGraphics
     * 
     */
    public mouseEvents(BranchGroup group,Canvas3D canvas,CreateGraphics parent)

{

    this.Cmain = parent;
    canvas3d = canvas;
    pickCanvas = new PickCanvas(canvas, group);

    pickCanvas.setMode(PickCanvas.GEOMETRY);
    
    
    canvas.addMouseListener(this);
    canvas.addMouseMotionListener(this);
}

    /**
     * 
     * Get coordinates of the mouse clicked on 3D canvas.
     * 
     * @param event
     * 
     * 
     * @return Coordinates
     */
    public Point3d getPosition(MouseEvent event) {
		Point3d eyePos = new Point3d();
		Point3d mousePos = new Point3d();
		canvas3d.getCenterEyeInImagePlate(eyePos);
		canvas3d.getPixelLocationInImagePlate(event.getX(),
                       event.getY(), mousePos);
		Transform3D transform = new Transform3D();
		canvas3d.getImagePlateToVworld(transform);
		transform.transform(eyePos);
		transform.transform(mousePos);
		Vector3d direction = new Vector3d(eyePos);
		direction.sub(mousePos);
		// three points on the plane
		Point3d p1 = new Point3d(.5, -.5, .5);
		Point3d p2 = new Point3d(.5, .5, .5);
		Point3d p3 = new Point3d(-.5, .5, .5);
		Transform3D currentTransform = new Transform3D();
//		box.getLocalToVworld(currentTransform);
		currentTransform.transform(p1);
		currentTransform.transform(p2);
		currentTransform.transform(p3);		
		Point3d intersection = getIntersection(eyePos, mousePos,
                        p1, p2, p3);
		currentTransform.invert();
		currentTransform.transform(intersection);
		return intersection;		
	}
    
    
    /**
     * @see mouseEvents#getPosition(java.awt.event.MouseEvent) 
     * 
     * @param line1
     * @param line2
     * @param plane1
     * @param plane2
     * @param plane3
     * 
     * @return coordinates
     */
Point3d getIntersection(Point3d line1, Point3d line2, 
			Point3d plane1, Point3d plane2, Point3d plane3) {
		Vector3d p1 = new Vector3d(plane1);
		Vector3d p2 = new Vector3d(plane2);
		Vector3d p3 = new Vector3d(plane3);
		Vector3d p2minusp1 = new Vector3d(p2);
		p2minusp1.sub(p1);
		Vector3d p3minusp1 = new Vector3d(p3);
		p3minusp1.sub(p1);
		Vector3d normal = new Vector3d();
		normal.cross(p2minusp1, p3minusp1);
		// The plane can be defined by p1, n + d = 0
		double d = -p1.dot(normal);
		Vector3d i1 = new Vector3d(line1);
		Vector3d direction = new Vector3d(line1);
		direction.sub(line2);
		double dot = direction.dot(normal);
		if (dot == 0) return null;
		double t = (-d - i1.dot(normal)) / (dot);
		Vector3d intersection = new Vector3d(line1);
		Vector3d scaledDirection = new Vector3d(direction);
		scaledDirection.scale(t);
		intersection.add(scaledDirection);
		Point3d intersectionPoint = new Point3d(intersection);
		return intersectionPoint;
	}

/**
 * On mouse click this event fires.
 * 
 * @param e 
 */
@Override
public void mouseClicked(MouseEvent e)

{
    pickCanvas.setShapeLocation(e);
    PickResult result = pickCanvas.pickClosest();

    if (result == null) {
        
        //JOptionPane.showMessageDialog(null,"In order to add atom click on the existing atom","Info",JOptionPane.INFORMATION_MESSAGE);
        this.Cmain.Mmain.log.append("!!! In order to add atom click on the existing atom...\n");

    } else {

       p = (Primitive)result.getNode(PickResult.PRIMITIVE);

       s = (Shape3D)result.getNode(PickResult.SHAPE3D);

       if (p != null) {
          
          this.Cmain.CreateSphere(p);
           
       } else if (s != null) {
           
       } else{

          System.out.println("null");

       }
    }
}
}
