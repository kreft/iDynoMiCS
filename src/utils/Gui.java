package utils;

import java.awt.BorderLayout;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;

import javax.swing.*;


@SuppressWarnings("serial")
public class Gui extends JPanel {
   private JTextArea textArea = new JTextArea(20, 60);
   private TextAreaOutputStream taOutputStream = new TextAreaOutputStream(
         textArea, " ");
   
   public static void openGui(String title) {
      try {
		SwingUtilities.invokeAndWait(new Runnable() {
		      public void run() {
		         Gui.createAndShowGui(title);
		      }
		   });
		} catch (InvocationTargetException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

   public Gui() {
      setLayout(new BorderLayout());
      add(new JScrollPane(textArea, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, 
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER));
      System.setOut(new PrintStream(taOutputStream));
   }

   public static void createAndShowGui(String title) {
      JFrame frame = new JFrame(title);
      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      frame.getContentPane().add(new Gui());
      frame.pack();
      frame.setLocationRelativeTo(null);
      frame.setVisible(true);
   }
}