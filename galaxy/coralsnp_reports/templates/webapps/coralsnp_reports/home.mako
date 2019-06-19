<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <a href="http://baumslab.org/" target="_blank"><img src="${h.url_for('/static/images/coral.jpg')}" width="500" height="344"></a>
        <br>
        <p>
        </p>
        <p>
Using molecular techniques to answer fundamental questions about evolution and ecology to guide coral reef conservation efforts.
        </p>
    </div>
</div>
