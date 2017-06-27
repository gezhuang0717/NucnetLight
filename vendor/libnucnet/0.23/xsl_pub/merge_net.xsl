<?xml version="1.0" encoding="iso-8859-1" ?>

<!--////////////////////////////////////////////////////////////////////////////
// <file type = "public">
//
//   <license>
//      See the README.txt file in this directory for copyright and license
//      information.
//   </license>
//
//   <description>
//     <abstract>
//       Stylesheet to merge an input Libnucnet__Nuc xml file and an 
//       input Libnucnet__Reac xml file to form an appropriate input
//       Libnucnet__Net xml file.
//     </abstract>
//     <keywords>
//       merge, nuc, reac, net, xml, file
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2007/08/06" />
//     </current>
//     <previous>
//     </previous>
//   </authors>
//
//   <compatibility>
//     Tested with xsltproc using libxml 20629, libxslt 10121 and
//     libexslt 813 xsltproc was compiled against libxml 20629,
//     libxslt 10121 and libexslt 813 libxslt 10121 was compiled against
//     libxml 20629 libexslt 813 was compiled against libxml 20629
//   </compatibility>
//
// </file>
/////////////////////////////////////////////////////////////////////////////-->

<xsl:stylesheet
  version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:nuc="http://libnucnet.sf.net/xsd_pub/2011-05-12/libnucnet__nuc/"
  xmlns:reac="http://libnucnet.sf.net/xsd_pub/2011-05-12/libnucnet__reac/"
>

<xsl:param name="reac_doc"/>

<xsl:template match="/">

  <xsl:element
     name="nuclear_network"
  >

    <xsl:element name="nuclear_data">

      <xsl:for-each select="child::nuclear_data">
        <xsl:apply-templates/>
      </xsl:for-each>

    </xsl:element>

    <xsl:element name="reaction_data">

      <xsl:for-each select="document($reac_doc)//child::reaction_data">
        <xsl:apply-templates/>
      </xsl:for-each>

    </xsl:element>

  </xsl:element>

</xsl:template>

<xsl:template match="node() | @* | comment() | processing-instruction()">
  <xsl:copy>
    <xsl:apply-templates select="@* | node()"/>
  </xsl:copy>
</xsl:template>

</xsl:stylesheet>

