package molecular_modeling;

import com.sun.j3d.utils.behaviors.vp.OrbitBehavior;
import com.sun.j3d.utils.universe.SimpleUniverse;
import com.sun.j3d.utils.geometry.Primitive;
import com.sun.j3d.utils.geometry.Sphere;
import com.sun.j3d.utils.image.TextureLoader;
import com.sun.j3d.utils.picking.PickTool;
import javax.media.j3d.*;
import javax.vecmath.Color3f;
import javax.vecmath.Point3f;
import java.awt.GraphicsConfiguration;
import java.awt.Toolkit;
import java.util.Enumeration;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3f;

/**
 *
 * @author Syed ahsan ali
 * 
 * This class is responsible for creating atoms and links between them on the 3D canvas.
 * 
 * Controls all the function regarding graphics.
 */
public class CreateGraphics {
    static BranchGroup scene;
    static Canvas3D canvas3D;
    static OrbitBehavior orbit;
    static TransformGroup[] trans = new TransformGroup[1000];
    static Shape3D[] dotShape = new Shape3D[1000];
    static Vector3f[] vertex = new Vector3f[1000];
    static String[] linesInfo = new String[1000];
    static int nextSphere = 0, HighlightedID = -1;
    static int nextLine = 0;
    static SimpleUniverse u;

    public static molecularModeling Mmain;
    static String Direction = null;
    
    TransformGroup createTG(float x, float y, float z) {
    Vector3f position = new Vector3f(x, y, z);
    Transform3D translate = new Transform3D();
    translate.set(position);
    TransformGroup TG = new TransformGroup(translate);
    return TG;
  }
    
 static TransformGroup createTg(float x, float y, float z) {
   Vector3f position = new Vector3f(x, y, z);
    Transform3D translate = new Transform3D();
    translate.set(position);
    TransformGroup TG = new TransformGroup(translate);
    return TG;
  }
 
  static Appearance createMatAppear(Color3f dColor, Color3f sColor, float shine) {

    Appearance appear = new Appearance();
    Material material = new Material();
    material.setDiffuseColor(dColor);
    material.setSpecularColor(sColor);
    material.setShininess(shine);
    appear.setMaterial(material);

    return appear;
  }

    /**
     * 
     * Calls only once after starting the program to create 3D canvas and initial atom and setting 
     * up the 3D environment, also register mouse events by calling MouseEvents class.
     * 
     * @see molecularModeling#StartSceneGraphWindow()
     *
     * @param parent molecularModeling class reference
     * 
     * @return 3D Canvas Object
     * 
     */
      public Canvas3D CreateGraphics(molecularModeling parent) {
        this.Mmain = parent;
        
        GraphicsConfiguration gc =SimpleUniverse.getPreferredConfiguration();
        canvas3D = new Canvas3D(gc);
        
        Toolkit.getDefaultToolkit().sync();
        scene = new BranchGroup();
        
        scene.setCapability(BranchGroup.ALLOW_DETACH);
        scene.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
        scene.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
        scene.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
        scene.setCapability(BranchGroup.ALLOW_BOUNDS_READ);
        scene.setCapability(BranchGroup.ALLOW_BOUNDS_WRITE);
        
        trans[nextSphere] = createTG(0.0f, 0.0f, 0.0f);
        vertex[nextSphere] = new Vector3f(0.0f,0.0f,0.0f);
        //trans.setTransform(createTG(0.0f, 0.0f, 0.0f));
        scene.addChild(trans[nextSphere]);
        Sphere sp = new Sphere(0.035f, Sphere.GENERATE_NORMALS, 60,
        createMatAppear(new Color3f(1.0f, 0.0f, 0.0f),new Color3f(1.0f, 1.0f, 1.0f), 50.0f));
        sp.setName("S-"+nextSphere);
        sp.setUserData("0");
        trans[nextSphere].setName("S-"+nextSphere);
        trans[nextSphere].addChild(sp);
        enablePicking(trans[nextSphere]);
        nextSphere++;

        mouseEvents mE=new mouseEvents(scene,canvas3D,this);

    AmbientLight lightA = new AmbientLight();
    lightA.setInfluencingBounds(new BoundingSphere());
    scene.addChild(lightA);

    DirectionalLight lightD1 = new DirectionalLight();
    lightD1.setInfluencingBounds(new BoundingSphere());
    Vector3f direction = new Vector3f(-1.0f, -1.0f, -1.0f);
    direction.normalize();
    lightD1.setDirection(direction);
    lightD1.setColor(new Color3f(1.0f, 1.0f, 1.0f));
    scene.addChild(lightD1);
    
    DirectionalLight lightD2 = new DirectionalLight();
    lightD2.setInfluencingBounds(new BoundingSphere());
    Vector3f direction2 = new Vector3f(1.0f, 1.0f, 1.0f);
    direction.normalize();
    lightD2.setDirection(direction2);
    lightD2.setColor(new Color3f(1.0f, 1.0f, 1.0f));
    scene.addChild(lightD2);
    
    Background background = new Background(new Color3f(0.752941f, 0.752941f, 0.752941f));

    background.setApplicationBounds(new BoundingSphere(new Point3d(0.0,0.0,0.0),Double.MAX_VALUE));
    scene.addChild(background);
    
     u = new SimpleUniverse(canvas3D);
    orbit = new OrbitBehavior(canvas3D, OrbitBehavior.REVERSE_ROTATE);
    orbit.setSchedulingBounds(new BoundingSphere());
    
    u.getViewingPlatform().setViewPlatformBehavior(orbit);
     u.getViewingPlatform().setNominalViewingTransform();
     u.getViewingPlatform().getViewPlatform().setActivationRadius(300);
    u.getViewer().getView().setBackClipDistance ( 300.0 );
    u.getViewer().getView().setLocalEyeLightingEnable(true);
    rEstScene();
    
    return canvas3D;
    
  }

    /**
     * 
     * This function checks whether the links between two atoms exists or not, 
     *  if exists return id, if not return -1.
     *
     * @param id1 Id of one atom
     * @param id2 Id of the other atom
     * 
     * @return int
     * 
     */
      public static int lineExistbetween(int id1,int id2)
  {
      String value = id1+"-"+id2;
      for(int i=0;i<linesInfo.length;i++)
        if(linesInfo[i] != null)
          if(linesInfo[i].equals(value))
              return i;
      
      value = id2+"-"+id1;
      
      for(int i=0;i<linesInfo.length;i++)
          if(linesInfo[i] != null)
                if(linesInfo[i].equals(value))
                    return i;
      
      return -1;
  }
  
  

    /**
     * 
     * Check whether the atom exists on specified location
     * where new item is going to be drawn 
     * returns id of the atom, if the atom exists 
     * and return -1 if there is no atom exists before.
     * 
     * @see CreateGraphics#CreateSphere(com.sun.j3d.utils.geometry.Primitive)
     * 
     * @param i_vertex Coordinates of an atom
     * 
     * @return int
     * 
     */
      public static int SphereExistAtLocation(Vector3f i_vertex)
  {
      for(int i=0;i<vertex.length;i++)
          if(vertex[i] != null)
            if(vertex[i].x == i_vertex.x && vertex[i].y == i_vertex.y && vertex[i].z == i_vertex.z)
                return i;
      
      
      return -1;
  }
 

    /**
     * 
     * Sets rotation of the molecule according to the given points or 
     * center of the mass while rotating with the mouse.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * @see molecularModeling#centerofmassbtnActionPerformed(java.awt.event.ActionEvent) 
     * 
     * @param x Position x-axis
     * @param y Position y-axis
     * @param z Position z-axis
     * 
     */
      public static void setCenterOfMass(float x,float y,float z)
  {
      Point3d point3D= new Point3d(x, y, z);
      orbit.setRotationCenter(point3D);  
  }
  
  

    /**
     * remove the very latest atom drawn on the canvas.
     */
      public static void undoScene()
  {
      
      nextSphere--;
      if(nextSphere > 0)
      {
          scene.detach();
          scene.removeChild(trans[nextSphere]);
          //int nextline = nextSphere - 1;
          nextLine--;
          scene.removeChild(dotShape[nextLine]);
          rEstScene();
      }
      
  }
  
  

    /**
     * 
     * This function sets the capability of the 3D elements on the canvas 
     * to detect exactly the element when user clicks on it, means intersection.
     * 
     * @see CreateGraphics#CreateSphere(com.sun.j3d.utils.geometry.Primitive)
     * @see CreateGraphics#RedrawSphere(int, float, float, float)
     * @see CreateGraphics#DrawHighlighted(int, java.lang.String[], float, javax.vecmath.Color3f)
     * 
     * @param node Object
     * 
     */
      public static void enablePicking(Node node) {

    node.setPickable(true);

    node.setCapability(Node.ENABLE_PICK_REPORTING);

    try {

       Group group = (Group) node;

       for (Enumeration e = group.getAllChildren(); e.hasMoreElements();) {

          enablePicking((Node)e.nextElement());

       }

    }

    catch(ClassCastException e) {

        // if not a group node, there are no children so ignore exception

    }

    try {

          Shape3D shape = (Shape3D) node;
          
          PickTool.setCapabilities(node, PickTool.INTERSECT_FULL);

          for (Enumeration e = shape.getAllGeometries(); e.hasMoreElements();) {

             Geometry g = (Geometry)e.nextElement();

             g.setCapability(g.ALLOW_INTERSECT);

          }

       }

    catch(ClassCastException e) {

       // not a Shape3D node ignore exception

    }

}
    /**
     * 
     * This function used to create atom on the 3D canvas by setting all
     * the required properties, according to the requirements.
     * 
     * @see mouseEvents#mouseClicked(java.awt.event.MouseEvent)
     * 
     * @param shape Clicked Object
     * 
     */
      public void CreateSphere(Primitive shape)
  {
     //JOptionPane.showMessageDialog(null, shape + " / look");
        Color3f white = new Color3f(1.0f, 1.0f, 1.0f);
        Color3f blue = new Color3f(0.0f, 0.0f, 1.0f);
        Color3f yellow = new Color3f(1.0f, 1.0f, 0.0f);
        Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
        Color3f black = new Color3f(0.0f, 0.0f, 0.0f);
        //nextSphere++;
        
          Transform3D td = new Transform3D();
          Vector3f position = new Vector3f();
          shape.getLocalToVworld(td);
          td.get(position);
         int nextCase = Integer.parseInt(shape.getUserData().toString());
         //JOptionPane.showMessageDialog(null, nextCase);
         while(nextCase < 6)
            {

            
         Vector3f NpoS = getXYZcorD(nextCase,position);

         int i_id2 = SphereExistAtLocation(NpoS);
         //System.out.println("iiiiiii222222 "+i_id2);
         if( i_id2 != -1)
         {
             String[] val = shape.getName().split("-");
             int i_id1 = Integer.parseInt(val[1].toString());
             //System.out.println("iiiiiii1111111 "+i_id1);
             nextCase++;
             shape.setUserData(String.valueOf(nextCase));
         }
         else
         {
          nextCase++;
          shape.setUserData(String.valueOf(nextCase));
          //String[] getOld = shape.getName().split("-");
        
        trans[nextSphere] = createTG(NpoS.x, NpoS.y, NpoS.z);
        vertex[nextSphere] = new Vector3f(NpoS.x, NpoS.y, NpoS.z);
        //trans.setTransform(createTG(NpoS.x, NpoS.y, NpoS.z));
        scene.detach();
        scene.addChild(trans[nextSphere]);
        Sphere sp = new Sphere(0.035f, Sphere.GENERATE_NORMALS, 60,
        createMatAppear(red, white, 50.0f));
        sp.setName("S-"+nextSphere);
        sp.setUserData("0");
        trans[nextSphere].setName("S-"+nextSphere);
        //trans[nextSphere].setUserData(Direction);
        trans[nextSphere].addChild(sp);
        enablePicking(trans[nextSphere]);
        
            int index = nextLine;
            Point3f CuRrPoint = new Point3f(NpoS.x, NpoS.y, NpoS.z);
            Point3f PrEvpoint = new Point3f(position.x,position.y,position.z);
            DrawLinE(PrEvpoint,CuRrPoint,index);
            rEstScene();
            String[] id = shape.getName().split("-");
            linesInfo[nextLine] = id[1].toString()+"-"+nextSphere;
            this.Mmain.addTableRows(shape.getName(), sp.getName(), position, NpoS);
            nextSphere++;
            nextLine++;
            
            break;
         }
          
         }
  }
  
  

    /**
     * 
     * Function used to redraw sphere if its deleted or its coordinates changes.
     * 
     * @see molecularModeling#ReadPMFile(java.lang.String)
     * 
     * @param Id Id of an atom
     * @param x Points on x-axis
     * @param y Points on y-axis
     * @param z Points on z-axis
     * 
     */
      public static void RedrawSphere(int Id,float x,float y,float z)
  {
        Color3f white = new Color3f(1.0f, 1.0f, 1.0f);
        Color3f blue = new Color3f(0.0f, 0.0f, 1.0f);
        Color3f yellow = new Color3f(1.0f, 1.0f, 0.0f);
        Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
        Color3f black = new Color3f(0.0f, 0.0f, 0.0f);
        
        trans[Id] = createTg(x, y, z);
        vertex[Id] = new Vector3f(x, y, z);
        //trans.setTransform(createTG(NpoS.x, NpoS.y, NpoS.z));
        scene.detach();
        scene.addChild(trans[Id]);
        Sphere sp = new Sphere(0.035f, Sphere.GENERATE_NORMALS, 60,
        createMatAppear(red, white, 50.0f));
        sp.setName("S-"+Id);
        //sp.setUserData(neigh);
        trans[Id].setName("S-"+Id);
        trans[Id].setUserData("null");
        trans[Id].addChild(sp);
        enablePicking(trans[Id]);
        nextSphere = Id + 1;
  }

    /**
     * To set the discription of the atom or sphere
     * Discription: information about, who its conected to? what its id etc
     * 
     * important!.
     *
     * @param id Id of an atom
     * @param Discription Information about an atom
     * 
     */
      public static void SetDiscription(int id,String Discription)
    {
      trans[id].getChild(0).setUserData(Discription);
  }
  
    /**
     * 
     * Saves the information about the connection b/w two atoms.
     *
     * @param ID Id of the connection
     * @param a  Id of one atom
     * @param b  Id of second atom
     * 
     */
    public static void SetLineArray(int ID,int a, int b)
  {
      linesInfo[ID] = a+"-"+b;
      nextLine = ID + 1;
  }
  
    /**
     * 
     * Assign the coordinates simultaneously while users clicks on the canvas 
     * in order to make the molecule first time.
     *
     * @param rAnGE Direction of the atom
     * @param poS Position in vector form
     * 
     * @return Vector3f Initial coordinates for the atom
     * 
     */
  private Vector3f getXYZcorD(int rAnGE,Vector3f poS)
  {
      Vector3f cord = poS; 
      switch(rAnGE)
              {
                case 0: cord = new Vector3f(poS.x + 0.15f,poS.y,poS.z);
                    break;
                case 1: cord = new Vector3f(poS.x - 0.15f,poS.y,poS.z);
                    break;
                case 2: cord = new Vector3f(poS.x,poS.y + 0.15f,poS.z);
                    break;
                case 3: cord = new Vector3f(poS.x ,poS.y - 0.15f,poS.z);
                    break;
                case 4: cord = new Vector3f(poS.x ,poS.y,poS.z - 0.15f);
                    break;
                case 5: cord = new Vector3f(poS.x ,poS.y,poS.z + 0.15f);
                    break;
              }
      return cord;
  }

    /**
     * 
     * This function draw a connection or link between two atoms 
     * based on coordinates of the two atoms.
     * 
     * @see CreateGraphics#CreateSphere(com.sun.j3d.utils.geometry.Primitive)
     * @see CreateGraphics#changeVertex(int, float, float, float, java.lang.String)
     * 
     * @param prEvPoint Coordinates of one atom
     * @param cuRrPoint Coordinates of second atom
     * @param nextline  Number of connectivity
     * 
     */
      public static void DrawLinE(Point3f prEvPoint,Point3f cuRrPoint, int nextline)
  {
      Point3f[] dotPts = new Point3f[2];
      Color3f black = new Color3f(0.0f, 0.0f, 0.0f);
        ColoringAttributes ca = new ColoringAttributes(black,
        ColoringAttributes.SHADE_GOURAUD);
        dotPts[0] = prEvPoint;
        dotPts[1] = cuRrPoint;
        LineArray dot = new LineArray(2, LineArray.COORDINATES);
        dot.setCoordinates(0, dotPts);
        LineAttributes dotLa = new LineAttributes();
        dotLa.setLineWidth(4.0f);
        dotLa.setLinePattern(LineAttributes.PATTERN_SOLID);
        Appearance dotApp = new Appearance();
        dotApp.setLineAttributes(dotLa);
        dotApp.setColoringAttributes(ca);
        dotShape[nextline] = new Shape3D(dot, dotApp);
        dotShape[nextline].setName("Sline-"+nextline);
        scene.addChild(dotShape[nextline]);
  }

    /**
     * Whenever any atom changes its position, this function is to be called to remove 
     * atom from its current position and draw it to the new location
     * Managing all the links it has and the information it has in the description.
     * 
     * @see molecularModeling#CalculateGeometryActionPerformed(java.awt.event.ActionEvent)
     * @see molecularModeling#addTablelistener()
     * 
     * @param Id Id of the atom to be moved
     * @param x New points on x-axis
     * @param y New points on y-axis
     * @param z New points on z-axis
     * @param link Existing links of an atom
     * 
     */
      public static void changeVertex (int Id,float x,float y,float z,String link)
  {
        Color3f white = new Color3f(1.0f, 1.0f, 1.0f);
        Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
        scene.detach();
        scene.removeChild(trans[Id]);

      String[] node = link.split(",");
      System.out.println(node[0]);
      
      String userdata =  trans[Id].getChild(0).getUserData().toString();
      
       trans[Id] = createTg(x, y, z);
       vertex[Id] = new Vector3f(x, y, z);

        scene.addChild(trans[Id]);
        Sphere sp = new Sphere(0.035f, Sphere.GENERATE_NORMALS, 60,
        createMatAppear(red, white, 50.0f));
        sp.setName("S-"+Id);
        sp.setUserData(userdata);
        trans[Id].setName("S-"+Id);
        trans[Id].addChild(sp);
        enablePicking(trans[Id]);
       
      for(int i=0;i < node.length; i++)
      {
          int inlink = Integer.parseInt(node[i].trim());
          int lineID = lineExistbetween(Id,inlink);
          if(lineID != -1)
          {
            scene.removeChild(dotShape[lineID]);
            Point3f CuRrPoint = new Point3f(x, y, z);
            Point3f PrEvpoint = new Point3f(vertex[inlink].x,vertex[inlink].y,vertex[inlink].z);
            DrawLinE(PrEvpoint,CuRrPoint,lineID);
          }
      }
      
      rEstScene();
      
  
  }

    /**
     * Changes the color of the atom when user selects any atom 
     * in the main table of the atoms.
     * 
     *@see molecularModeling#addTablelistener()
     * 
     *@param ID Id of an atom to be highlighted
     * 
     * 
     */
      public static void HighlightSphere(int ID)
  {
      
      Color3f black = new Color3f(0.0f, 0.0f, 0.0f);
      Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
      String[] data = new String[3];
      
      if(HighlightedID == -1)
      {
          HighlightedID = ID;
          
          data[0] = trans[ID].getChild(0).getName();
          data[1] = trans[ID].getChild(0).getUserData().toString();
          data[2] = trans[ID].getName();
          
          DrawHighlighted(ID,data,0.040f,black);
      }
      else
      {
          data[0] = trans[HighlightedID].getChild(0).getName();
          data[1] = trans[HighlightedID].getChild(0).getUserData().toString();
          data[2] = trans[HighlightedID].getName();
          
          DrawHighlighted(HighlightedID,data,0.035f,red);
          
          data[0] = trans[ID].getChild(0).getName();
          data[1] = trans[ID].getChild(0).getUserData().toString();
          data[2] = trans[ID].getName();
          
          DrawHighlighted(ID,data,0.040f,black);
          
          HighlightedID = ID;
      }
  }
  
    /**
     * 
     *@see CreateGraphics#HighlightSphere(int)
     * 
     */
  private static void DrawHighlighted(int Id,String[] InfoData,float size,Color3f color)
  {
      Color3f white = new Color3f(1.0f, 1.0f, 1.0f);
      
      scene.detach();
     scene.removeChild(trans[Id]);
    
       trans[Id] = createTg(vertex[Id].x, vertex[Id].y, vertex[Id].z);
    
        scene.addChild(trans[Id]);
        Sphere sp = new Sphere(size, Sphere.GENERATE_NORMALS, 60,
        createMatAppear(color, white, 50.0f));
        sp.setName(InfoData[0]);
        sp.setUserData(InfoData[1]);
        trans[Id].setName(InfoData[2]);
       // trans[Id].setUserData(InfoData[3]);
        trans[Id].addChild(sp);
        enablePicking(trans[Id]);
        rEstScene();
  }

    /**
     * rest of the setting of the scene when any change occure on the 3D canvas.
     */
      public static void rEstScene()
  {
    u.addBranchGraph(scene);
    canvas3D.getView().repaint();
  }
      
    /**
     * reset the molecule to its original position.
     */
      public static void ResetScene()
  {
      u.getViewingPlatform().setNominalViewingTransform();
  }

    /**
     * cleans the 3D canvas or scene graph before reading an 
     * existing structure from the file
     * reset every needed settings to completely start a new one.
     */
      public static void CleanSceneGraph()
  {
      scene.detach();
      for(int i=0;i<trans.length;i++)
      {
          scene.removeChild(trans[i]);
          trans[i] = null;
          vertex[i] = null;
      }
      for(int j=0;j<dotShape.length;j++)
      {
          scene.removeChild(dotShape[j]);
          dotShape[j] =  null;
      }
      rEstScene();
      for(int k=0;k<linesInfo.length;k++)
      {
          linesInfo[k] = null;
      }
      nextSphere = 0;
      HighlightedID = -1;
      nextLine = 0;
      ResetScene();
  }
}