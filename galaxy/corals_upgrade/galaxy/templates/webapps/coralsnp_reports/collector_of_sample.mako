<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Collector of sample with affy id ${affy_id}
        </h3>
        <table align="center" class="colored">
            %if len(collectors) == 0:
                <tr>
                    <td colspan="2">
                        This sample has no collector
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Last Name</td>
                    <td>First Name</td>
                    <td>Organization</td>
                    <td>Email</td>
                </tr>
                <% ctr = 0 %>
                %for collector in collectors:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${collector[0]}</td>
                        <td>${collector[1]}</td>
                        <td>${collector[2]}</td>
                        <td>${collector[3]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

